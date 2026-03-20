# xlsx_compiler.py

"""
Name:       xlsx_compiler.py
Author:     CAG
Version:    1.5
Date:       2026/03/20
"""

#%% Imports

import utils.qc_ops as qc
import pandas as pd
import numpy as np
import logging
from typing import List

logger = logging.getLogger(__name__)

#%% Classes

class set_compiler:
    """
    set_compiler class: Combines parsed files into a single DataFrame and
    runs all transformation steps, ending with above-cutoff filtering
    and writing outputs.
    """

    def __init__(
        self, 
        parsed_files: List[pd.DataFrame],
        settings: dict,
        base_df: pd.DataFrame | None = None # optional base dataframe
        ) -> None:
        self.base_df = base_df
        self.df: pd.DataFrame = self._combine_new(parsed_files)
        self._sort()
        self._base_qc()
        self._optional_qc(settings)
        self._above_cutoff_df() # after ldist_check so it contains ldist_rerun


    #####################
    # init methods
    #####################

    def _combine_new(self, parsed_files: List[pd.DataFrame]) -> pd.DataFrame:
        """
        Combine parsed DataFrames into a single DataFrame.
        """
        if not parsed_files:
            return pd.DataFrame()

        df_new = pd.concat(parsed_files, ignore_index=True)

        # no base → original behavior
        if self.base_df is None or self.base_df.empty:
            return df_new

        # base present → apply override logic
        return self._combine_with_base(self.base_df, df_new)


    def _combine_with_base(self, df_base: pd.DataFrame, df_new: pd.DataFrame) -> pd.DataFrame:
        """
        Compare base csv data to new analyses set.
        Keep all new analysis samples, and any non-overlapping base samples.

        Assumption: any re-analyzed (or newly modified) run is the "correct" one.
        Identical sample/run in new analysis OVERWRITES sample/run in base.
        """
        # These columns capture a single sample from a given run
        id_cols = ["samp_group", "samp_name", "samp_date", "run_number"]
        
        # Standardize data types before comparison
        for df in [df_base, df_new]:
            # ID columns
            df["samp_group"] = df["samp_group"].astype(str)
            df["samp_name"] = df["samp_name"].astype(str)
            df["samp_date"] = df["samp_date"].apply(qc.normalize_date)
            df["run_number"] = df["run_number"].astype(str)
            
            # Other pivot/grouping columns downstream...
            df["filename"] = df["filename"].astype(str)
            df["bc_name"] = df["bc_name"].astype(str)
            df["bc_seq"] = df["bc_seq"].astype(str)
        
        # ID discrete samples... 
        ids_new = pd.MultiIndex.from_frame(df_new[id_cols].drop_duplicates())
        base_row_ids = pd.MultiIndex.from_frame(df_base[id_cols])
        
        # Mask - base data samples not in new data
        base_keep_mask = ~base_row_ids.isin(ids_new)

        # Combine base-unique and total new samples to single df 
        return pd.concat([df_base.loc[base_keep_mask], df_new], ignore_index=True)


    def _sort(self) -> None:
        """
        Rearrange columns to prioritize certain metadata fields
        and sort by sample group, name, and date.
        """
        df = self.df.copy()
        cols = ['samp_group', 'samp_name', 'samp_date']

        # Move priority columns to the front
        df = df[[c for c in cols if c in df.columns] +
                [c for c in df.columns if c not in cols]]

        # Sort by ascending names and dates
        df = df.sort_values(
            by=cols,
            ascending=[True] * len(cols),
            kind='stable'
            ).reset_index(drop=True)

        self.df = df


    def _base_qc(self) -> None:
        """
        Run a series of basic qc operations
        """
        df = self.df.copy()

        # Make sure these are numeric
        cols = ["bc_count", "proportion", "input"]
        df[cols] = df[cols].apply(pd.to_numeric, errors='coerce').fillna(0)

        # Add ranks to 'unique' names
        df = qc.rank_uniques(df=df)

        # Fix instances of "293" from analyses using old ref. There is no 293. 
        df["bc_name"] = df["bc_name"].str.replace("SIVmac293M2", "SIVmac239M2", regex=False)

        # Flag short barcodes
        df = qc.flag_shorts(df = df)

        # Standardize hdist format across runs
        df["ldist_samp_lvl"] = df["ldist_samp_lvl"].str.replace(r'^(1-)|(-1)$', '', regex=True)
        df["ldist_samp_lvl"] = df["ldist_samp_lvl"].str.replace(r'^(Unique)', '', regex=True)
        df["ldist_samp_lvl"] = df["ldist_samp_lvl"].replace("NA-NA", None)

        # Flag same barcode across multiple samples in run (indicates shared, if same samp_group, may be expected)
        df["multi_idx"] = (
            df.groupby(["run_number", "samp_group", "bc_name"])["samp_name"]
            .transform(lambda s: s.nunique() > 1)
            )

        self.df = df


    def _optional_qc(self, settings):
        """
        Flesh this out to include
        - [X] Per-sample ldist checks (replace original hdist analysis)
        - [X] Per-group ldist checks
        - [ ] Comprehensive index-hopping/contamination assessment?
        """

        if settings["ldist_samp_lvl"] == True:
            self._ldist_qc(
                groupby_cols=['samp_group', 'samp_name', 'samp_date'],
                output_col='ldist_samp_lvl',
                drop_existing=True,
                settings=settings,
                )

        if settings["ldist_group_lvl"]:
            self._ldist_qc(
                groupby_cols=['samp_group'],
                aggr=True,
                output_col='ldist_group_lvl',
                settings=settings,
                )

        if settings["contam_check"]:
            self._ldist_qc(
                groupby_cols=None,
                aggr=True,
                output_col='ldist_all',
                settings=settings,
                )
            self.contam_check=True
        else:
            self.contam_check=False # to pass to io/write_sc_output


    def _ldist_qc(self, settings, groupby_cols, aggr=False, output_col='putative_parent', drop_existing=False):
        """
        Run ldist QC checks with flexible grouping and column assignments.

        Parameters:
        - groupby_cols (list): Columns to group the dataframe by.
        - aggr (bool): Whether to aggregate results inside flag_putative_parents.
        - output_col (str): Name of the output column to hold QC flags.
        - drop_existing (bool): Drop existing output_col column before recalculation.
        """
        if drop_existing and output_col in self.df.columns:
            self.df = self.df.drop(columns=[output_col])

        if groupby_cols:
            # run at given groupby level
            self.df = (
                self.df
                .groupby(groupby_cols, as_index=False)
                .apply(lambda g: qc.flag_putative_parents(
                    set_name=g.name,
                    df=g,
                    settings=settings,
                    aggr=aggr
                )).rename(columns={'putative_parent': output_col})
                )
        else:
            # run for entire set
            self.df = (qc.flag_putative_parents(
                set_name="Full Series",
                df = self.df,
                settings=settings,
                aggr=aggr
                ).rename(columns={'putative_parent': output_col})
                )


    def _above_cutoff_df(self) -> None:
        """
        Generate an above-cutoff DataFrame and recalc proportions
        using *post-cutoff* totals per (samp_name, samp_date, run_number, filename).
        """
        df_ac = self.df.loc[self.df["above_cutoff"]].copy()

        # Define sample identity
        sample_keys = ["samp_group","samp_name", "samp_date", "run_number", "filename"]

        # Get post-cutoff totals, including groups with NaNs
        df_ac["__denom"] = (
            df_ac.groupby(sample_keys, dropna=False)["bc_count"]
                .transform("sum")
                .astype(float)
        )

        # Safe division; make the denominator inspectable
        with np.errstate(divide="ignore", invalid="ignore"):
            df_ac["proportion"] = (df_ac["bc_count"] / df_ac["__denom"]).astype(float)

        df_ac["proportion"] = df_ac["proportion"].replace([np.inf, -np.inf], np.nan).fillna(0.0)

        # drop this bad boy
        df_ac = df_ac.drop(columns="__denom")

        self.df_ac = df_ac
    #################
    # external methods
    #################

    def group_mat_dict(self) -> None:
        """
        Format a dict of per-group above-cutoff sample matrices
        """
        df = self.df_ac.copy()
        df_dict = {k: v for k, v in df.groupby("samp_group")}
        mat_dict = {}

        for samp_group, samp_df in df_dict.items():
            
            # Check for duplicates in pivot keys
            key_cols = [
                "bc_name",
                "bc_seq",
                "samp_name",
                "samp_date",
                "run_number",
                "filename",
                "input",
            ]

            # Find duplicates
            dup_mask = samp_df.duplicated(subset=key_cols, keep=False)

            if dup_mask.any():
                print("\n===============================================")
                print(f"ERROR: Duplicate pivot key combinations in samp_group = {samp_group}")
                print("Offending rows:")
                print(samp_df.loc[dup_mask].sort_values(key_cols)[key_cols + ["proportion"]])
                print("===============================================")

                raise SystemExit(
                    f"Duplicate rows detected in samp_group={samp_group}; cannot pivot."
                )

            # Pivot into matrix
            df_mat = (
                samp_df.pivot(
                    index=["bc_name", "bc_seq"],
                    columns=["samp_name", "samp_date", "run_number", "filename", "input"],
                    values="proportion",
                )
                .fillna(0)
                .sort_index(axis=1, level=["samp_name", "samp_date"])
                .reset_index()
                )

            # Sort by all proportion columns
            cols_to_sort = df_mat.columns[2:]
            df_mat = df_mat.sort_values(
                by=list(cols_to_sort),
                ascending=[False] * len(cols_to_sort),
                na_position="last",
                ).reset_index(drop=True)

            mat_dict[samp_group] = df_mat
        
        return(df_dict, mat_dict)


    def contam_report(self) -> None:
        """
        Format a dict of per-group above-cutoff sample matrices
        """
        df = self.df_ac.copy()

        # Add new per-barcode columns
        df['n_groups'] = df.groupby('bc_name')['samp_group'].transform('nunique')
        df['n_samps'] = df.groupby('bc_name')['samp_name'].transform('nunique')

        # Add per-row columns
        df["input"] = df["input"].replace([np.inf, -np.inf, np.nan], 1)
        df["adj_count"] = df["proportion"] * df["input"]

        # --- Full matrix ---
        # Pivot to multi-index
        full_mat = (
            df.pivot(
                index=["bc_name", "bc_seq", "ldist_all", "n_groups", "n_samps"],
                columns=["samp_group", "samp_name", "samp_date", "run_number", "filename", "input"],
                values="proportion",
            )
            .fillna(0)
            .sort_index(axis=1, level=["samp_group","samp_name", "samp_date"])
            .reset_index()
            )

        # Sort by all proportion columns
        cols_to_sort = full_mat.columns[5:]
        full_mat = full_mat.sort_values(
            by=list(cols_to_sort),
            ascending=[False] * len(cols_to_sort),
            na_position="last",
            ).reset_index(drop=True)

        # --- Contamination report ---
        agg_df = (
            df.groupby(['bc_name', 'bc_seq', 'n_groups', 'ldist_all', 'short_bc', 'samp_group'], dropna=False)
            .agg(
                n_samp_per_group=('samp_name', 'nunique'),
                group_adj_count=('adj_count', 'sum')
            )
            .reset_index()
            )

        contam_df = (
            agg_df
            # Calculate total adjusted counts per animal by mapping sum back
            .assign(total_adj_count_animal=agg_df.groupby('samp_group')['group_adj_count'].transform('sum'))
            # Calculate proportion at barcode level relative to animal total
            .assign(group_prop=lambda d: d['group_adj_count'] / d['total_adj_count_animal'])
            # Filter only barcodes with more than one group
            .query('n_groups > 1')
            # Pivot with multiple value columns
            .pivot(
                index=["bc_name", "bc_seq", 'n_groups', 'ldist_all', 'short_bc'],
                columns="samp_group",
                values=["n_samp_per_group", "group_adj_count", "group_prop"]
            )
            .reset_index()
            # Sort by descending n_groups
            .sort_values(by='n_groups', ascending=False)
            )

        return(full_mat, contam_df)

#%% Versions
"""
v1.5
- Refactored qc_ops to pass settings dict instead of indiv params
v1.4
- Can now add non-redundant samples from a provided base csv
v1.3
- Passing params via settings dict, and accomodating flag_putative_parents modifications
v1.2 
- moved shared methods w/ parse_countdata to utils.qc_ops and output method to utils.io
- added optional qc methods ldist_samp_lvl, ldist_group_lvl, contam_check
v1.1 
- added functionality: modify "SIVmac293" >> "SIVmac239" in bc_names, to ensure accurate matrix rows w/ old studies using wrong name.
"""
