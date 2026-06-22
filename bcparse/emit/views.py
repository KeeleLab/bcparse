# views.py
"""
Name:      views.py
Author:     CAG
Version:    1.0
Date:       2026/03/26
Refactored: 2026/05/26

Functions for reshaping emitted dataframes into viewer-friendly formats.
Lifted from deprecated io.py.
"""

#%% Imports
from __future__ import annotations

import numpy as np
import pandas as pd

import bcparse.qc as qc
from bcparse.containers.runseries import RunSeries
from bcparse.containers.seqrun import SeqRun

#%%

def _legacy_bool_string(value) -> str:
    if pd.isna(value):
        return ""
    return "Yes" if qc.coerce_boolish(value) else ""


def _display_bool_or_blank(value):
    if pd.isna(value):
        return ""
    return qc.coerce_boolish(value)


def _format_legacy_date(value) -> str:
    date_val = pd.to_datetime(value, errors="coerce")
    if pd.isna(date_val):
        return ""
    return date_val.strftime("%Y-%m-%d")


def sidelong_tables(df: pd.DataFrame, legacy_format: bool = False):
    """
    Format viewer-oriented sidelong sample tables from a SeqRun-style dataframe.
    """
    if df.empty:
        return {}, {}

    group_key = "samp_group"
    split_key = "sample_id"
    header_keys = [
        "run_number",
        "run_name",
        "run_date",
        "samp_group",
        "samp_name",
        "samp_date",
        "input",
    ]

    # Keep sample blocks in a predictable viewer order within each group.
    df = df.sort_values(
        by=["samp_group", "samp_name", "samp_date"],
        ascending=[True, True, True],
        kind="stable",
    ).reset_index(drop=True)

    # Get samp_group-to-sample_id mapping
    groups = df.groupby(group_key)[split_key].unique().to_dict()

    # Split long df to dict with keys = sample_id
    split_dfs = {key: value.copy() for key, value in df.groupby(split_key, sort=False)}

    for split_id, split_df in split_dfs.items():
        # Prefer the native analysis-pipe parent column, but fall back to the
        # workbook-imported ldist column when compiling from Analysis.xlsx.
        if "putative_parent" in split_df.columns and split_df["putative_parent"].notna().any():
            parent_col = "putative_parent"
        elif "ldist_samp_lvl" in split_df.columns and split_df["ldist_samp_lvl"].notna().any():
            parent_col = "ldist_samp_lvl"
        else:
            parent_col = None

        #########
        # body df
        #########

        df_body = split_df[
            [
                "bc_name",
                "bc_seq",
                "bc_count",
                "proportion",
                "x_idx",
                "short_bc",
                "above_cutoff",
            ]
        ].copy()
        df_body["putative_parent"] = split_df[parent_col] if parent_col else pd.NA
        df_body = df_body[
            [
                "bc_name",
                "bc_seq",
                "bc_count",
                "proportion",
                "putative_parent",
                "x_idx",
                "short_bc",
                "above_cutoff",
            ]
        ]

        # Set conditions for grouping barcodes in sidelong table:
        c1 = df_body["above_cutoff"]
        c2 = df_body["putative_parent"].notna()
        c3 = df_body["bc_name"].str.startswith("Unique", na=False)
        c4 = df_body["bc_name"].str.startswith("Spike", na=False)

        # Bitpack
        df_body["flag"] = (
            c1.astype(int) * 8    # 8: above_cutoff
            + c2.astype(int) * 4  # 4: has parent
            + c3.astype(int) * 2  # 2: Unique
            + c4.astype(int)      # 1: Spike
        )

        df_body["rowgroup"] = df_body["flag"].map({
            8: 0,    # Above cutoff, no parent, named
            10: 0,   # Above cutoff, no parent, Unique
            9: 1,    # Above cutoff, no parent, Spike
            12: 1,   # Above cutoff, has parent, named
            13: 1,   # Above cutoff, has parent, Spike
            14: 1,   # Above cutoff, has parent, Unique
        }).fillna(2)

        df_body = df_body.sort_values(by=["rowgroup", "proportion"], ascending=[True, False]).reset_index(drop=True)

        # Recalculate group 0's proportions to total 1.
        ac_tot_prop = df_body.loc[df_body["rowgroup"] == 0, "proportion"].sum()
        if ac_tot_prop != 0:
            df_body.loc[df_body["rowgroup"] == 0, "proportion"] /= ac_tot_prop
        else:
            df_body.loc[df_body["rowgroup"] == 0, "proportion"] = 0

        changes = df_body.index[df_body["rowgroup"] != df_body["rowgroup"].shift()].tolist()
        if changes:
            changes.pop(0)

        pad = pd.DataFrame([[np.nan] * len(df_body.columns)], columns=df_body.columns)
        for index in sorted(changes, reverse=True):
            df_body = pd.concat([df_body.iloc[:index], pad, df_body.iloc[index:]], ignore_index=True)

        df_body = df_body.drop(["flag", "rowgroup", "above_cutoff"], axis=1)

        if legacy_format:
            df_body["x_idx"] = df_body["x_idx"].map(
                lambda value: "" if pd.isna(value) else f"{value:.3g}"
            )
            df_body["short_bc"] = df_body["short_bc"].astype(object)
            nan_mask = ~df_body["bc_name"].isna()

            df_body.loc[nan_mask, "short_bc"] = df_body.loc[
                nan_mask, "short_bc"
            ].map(_legacy_bool_string)

            df_body.loc[~nan_mask, "x_idx"] = ""
            df_body.loc[~nan_mask, "short_bc"] = ""

            df_body.loc[nan_mask, "putative_parent"] = (
                df_body.loc[nan_mask, "putative_parent"].fillna("NA-NA").map(lambda x: f"{x}" if x != "NA-NA" else x)
            )

            df_body = df_body.rename(
                columns={
                    "bc_name": "Barcode",
                    "bc_seq": "Sequence",
                    "bc_count": "Counts",
                    "proportion": "Proportion",
                    "putative_parent": "Hamming Dist to Another Sequence",
                    "x_idx": "Multi-group share",
                    "short_bc": "Short barcode?",
                }
            )
        else:
            df_body["short_bc"] = df_body["short_bc"].map(_display_bool_or_blank)

        new_cols = [f"{i}" for i in range(len(df_body.columns) + 1)]
        df_body = pd.concat([pd.DataFrame([df_body.columns], columns=df_body.columns), df_body]).reset_index(drop=True)

        for i in range(len(new_cols) - len(df_body.columns)):
            df_body[i] = np.nan

        df_body.columns = new_cols

        ###########
        # header df
        ###########

        header_dict = {}
        for key in header_keys:
            values = split_df[key].drop_duplicates()
            if len(values) != 1:
                raise ValueError(f"{key} has multiple values in {split_key}={split_id}: {values.tolist()}")
            header_dict[key] = values.iloc[0]

        for key in ["run_date", "samp_date"]:
            if key in header_dict:
                header_dict[key] = _format_legacy_date(header_dict[key])

        if legacy_format:
            header_data = [
                [
                    header_dict["samp_group"],
                    header_dict["samp_name"],
                    header_dict["samp_date"],
                    header_dict["samp_group"],
                    header_dict["samp_name"],
                    header_dict["samp_date"],
                ]
                + [np.nan] * (len(new_cols) - 6),
                [
                    "Sequencer Input",
                    "Total # Sequences",
                    "# of Barcodes (Excluding Uniques)",
                    "Proportion Uniques",
                ]
                + [np.nan] * (len(new_cols) - 4),
                [
                    header_dict["input"],
                    split_df["bc_count"].sum(),
                    split_df[~split_df["bc_name"].str.contains("Unique")].shape[0],
                    split_df[split_df["bc_name"].str.contains("Unique")]["bc_count"].sum() / split_df["bc_count"].sum(),
                ]
                + [np.nan] * (len(new_cols) - 4),
                [np.nan] * len(new_cols),
            ]
            df_head = pd.DataFrame(header_data, columns=new_cols)
        else:
            df_head = pd.DataFrame(list(header_dict.items()), columns=new_cols[:2])
            for col in new_cols[2:]:
                df_head[col] = np.nan

            pad = pd.DataFrame([[np.nan] * len(new_cols)], columns=new_cols)
            df_head = pd.concat([df_head, pad], ignore_index=True)

        split_dfs[split_id] = pd.concat([df_head, df_body], ignore_index=True)

    # Stitch per-sample sidelong tables back together by sample group.
    sidelongs = {}
    for group_name, split_ids in groups.items():
        dfs = [split_dfs[split_id] for split_id in split_ids if split_id in split_dfs]
        sidelongs[group_name] = pd.concat(dfs, axis=1) if dfs else pd.DataFrame()

    return sidelongs, split_dfs


def samp_group_matrices(seq_run: SeqRun, filt_ac=False, collapse_to_parent=False):
    """
    Return a per-samp_group dict of barcode proportion matrices,
    optionally filtered for above_cutoff.
    """
    if collapse_to_parent:
        df = pd.concat(
            [sample.collapse_to_parent().df for sample in seq_run.sample_list],
            ignore_index=True,
        )
    else:
        df = seq_run.df.copy()

    df[["run_date", "samp_date"]] = df[["run_date", "samp_date"]].apply(
        lambda col: col.apply(qc.normalize_date)
    )

    # Get distinct {samp_group: [idx_name(s)]}
    groups = df.groupby("samp_group")["idx_name"].unique().to_dict()

    # Initialize dict to hold per-samp_group results
    df_mats = {}

    for samp_group, idx_names in groups.items():
        # Optional cutoff filter
        if filt_ac:
            # Pull df copy
            df_filt = df[(df["idx_name"].isin(idx_names))].copy()

            # Conditions for grouping barcodes to match sidelong table logic:
            c1 = df_filt["above_cutoff"]
            c2 = df_filt["putative_parent"].notna()
            c3 = df_filt["bc_name"].str.startswith("Unique", na=False)
            c4 = df_filt["bc_name"].str.startswith("Spike", na=False)

            # Bitpack the flag columns using the same explicit scheme as sidelong_tables().
            df_filt["flag"] = (
                c1.astype(int) * 8
                + c2.astype(int) * 4
                + c3.astype(int) * 2
                + c4.astype(int)
            )

            # filter above cutoff rows
            df_filt = df_filt[df_filt["flag"].isin([8, 10])].copy() # Only include above_cutoff barcodes without parents, excluding Spikes but including Uniques, to match sidelong table logic.

            # recalculate proportions based on above_cutoff data
            df_filt["proportion"] = (
                df_filt.groupby("idx_name")["proportion"].transform(lambda x: x / x.sum())
            )
        else:
            df_filt = df[(df["idx_name"].isin(idx_names))]

        # Generate bc proportion matrix
        df_mat = (
            df_filt.pivot(
                index=["bc_name", "bc_seq"],
                columns=["samp_name", "samp_date", "idx_name", "input"],
                values="proportion",
            )
            .fillna(0)
            .sort_index(axis=1, level=["samp_name", "samp_date"])
            .reset_index()
        )

        cols_to_sort = df_mat.columns[2:]
        df_mat = df_mat.sort_values(
            by=list(cols_to_sort),
            ascending=[False] * len(cols_to_sort),
            na_position="last",
        ).reset_index(drop=True)

        df_mats[samp_group] = df_mat

    return df_mats


def group_mat_dict(df_ac: pd.DataFrame):
    """
    Format a dict of per-group above-cutoff sample matrices.
    """
    df = df_ac.copy()
    df[["run_date", "samp_date"]] = df[["run_date", "samp_date"]].apply(
        lambda col: col.apply(qc.normalize_date)
    )

    df_dict = {k: v for k, v in df.groupby("samp_group")}
    mat_dict = {}

    for samp_group, samp_df in df_dict.items():
        key_cols = [
            "bc_name",
            "bc_seq",
            "samp_name",
            "samp_date",
            "run_number",
            "filename",
            "input",
        ]

        dup_mask = samp_df.duplicated(subset=key_cols, keep=False)

        if dup_mask.any():
            dup_rows = samp_df.loc[dup_mask].sort_values(key_cols)[key_cols + ["proportion"]]
            raise ValueError(
                f"Duplicate rows detected in samp_group={samp_group}; cannot pivot.\n"
                f"Offending rows:\n{dup_rows.to_string(index=False)}"
            )

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

        cols_to_sort = df_mat.columns[2:]
        df_mat = df_mat.sort_values(
            by=list(cols_to_sort),
            ascending=[False] * len(cols_to_sort),
            na_position="last",
        ).reset_index(drop=True)

        mat_dict[samp_group] = df_mat

    return df_dict, mat_dict


def contam_report(df_ac: pd.DataFrame):
    """
    Format full-matrix and contamination report outputs for a compiled series.
    """
    df = df_ac.copy()
    df[["run_date", "samp_date"]] = df[["run_date", "samp_date"]].apply(
        lambda col: col.apply(qc.normalize_date)
    )

    df["n_groups"] = df.groupby("bc_name")["samp_group"].transform("nunique")
    df["n_samps"] = df.groupby("bc_name")["samp_name"].transform("nunique")

    df["input"] = df["input"].replace([np.inf, -np.inf, np.nan], 1)
    df["adj_count"] = df["proportion"] * df["input"]

    full_mat = (
        df.pivot(
            index=["bc_name", "bc_seq", "ldist_all", "n_groups", "n_samps"],
            columns=["samp_group", "samp_name", "samp_date", "run_number", "filename", "input"],
            values="proportion",
        )
        .fillna(0)
        .sort_index(axis=1, level=["samp_group", "samp_name", "samp_date"])
        .reset_index()
    )

    cols_to_sort = full_mat.columns[5:]
    full_mat = full_mat.sort_values(
        by=list(cols_to_sort),
        ascending=[False] * len(cols_to_sort),
        na_position="last",
    ).reset_index(drop=True)

    agg_df = (
        df.groupby(["bc_name", "bc_seq", "n_groups", "ldist_all", "short_bc", "samp_group"], dropna=False)
        .agg(n_samp_per_group=("samp_name", "nunique"), group_adj_count=("adj_count", "sum"))
        .reset_index()
    )

    contam_df = (
        agg_df.assign(total_adj_count_animal=agg_df.groupby("samp_group")["group_adj_count"].transform("sum"))
        .assign(group_prop=lambda d: d["group_adj_count"] / d["total_adj_count_animal"])
        .query("n_groups > 1")
        .pivot(
            index=["bc_name", "bc_seq", "n_groups", "ldist_all", "short_bc"],
            columns="samp_group",
            values=["n_samp_per_group", "group_adj_count", "group_prop"],
        )
        .reset_index()
        .sort_values(by="n_groups", ascending=False)
    )

    return full_mat, contam_df


def concat_series_runinfo(run_series: RunSeries) -> pd.DataFrame:
    """
    Concatenate per-run raw runinfo tables for compile output.
    """
    runinfo_dfs = []

    for seq_run in run_series.run_list:
        runinfo = seq_run.runinfo
        if runinfo is None:
            continue

        runinfo_df = runinfo.copy_raw_df()
        runinfo_df["Date"] = runinfo_df["Date"].apply(qc.normalize_date)

        if runinfo_df.empty:
            continue

        if runinfo_df.columns.size > 0 and str(runinfo_df.columns[0]).strip() == "Run Number":
            runinfo_df = runinfo_df.iloc[:, 1:].copy()

        runinfo_df = runinfo_df.mask(
            runinfo_df.map(lambda x: isinstance(x, str) and x.strip() == ""),
            np.nan,
        ).dropna(how="all")
        if runinfo_df.empty:
            continue

        run_number = seq_run.run_number
        if pd.isna(run_number):
            run_number = runinfo.run_number

        runinfo_df.insert(0, "run_number", run_number)
        runinfo_df.insert(0, "run_name", seq_run.run_name or runinfo.run_name)
        runinfo_dfs.append(runinfo_df)

    if not runinfo_dfs:
        return pd.DataFrame()

    return pd.concat(runinfo_dfs, ignore_index=True)

#%% Versions
"""
2026-05-26 - Refactored into bcparse package structure
"""
