# parse_countdata.py

"""
Name:       parse_countdata.py
Author:     CAG
Version:    1.9
Date:       2026/2/02
"""

# %% Imports

from collections import defaultdict

import numpy as np
import pandas as pd
import utils.qc_ops as qc
from rapidfuzz.distance import Levenshtein

# %% Class definition


class parse_countdata(object):
    """
    Class to process the counts dictionary {idx_seq:{bc_seq:count}} into class object with attributes:
        - rawdict: {idx_seq:{bc_seq:count}}, raw counts
        - bc_refdict_kseq: {bc_seq:bc_name}, for mapping known barcodes
        - idx_expect_kseq: {idx_seq:idx_name}, for mapping for expected indexes
        - idx_known: {idx_name:{bc_seq:count}}, for named indexes
        - idx_other: {idx_seq:{bc_seq:count}}, for unknown indexes
        - bc_set: df with non-redundant barcode information
        - df: final processed dataframe

    output methods:
    parse_countdata.samp_group_matrices()
    parse_countdata.samp_group_sidelong_tables()
    """

    def __init__(
        self,
        countdict: dict,
        runinfo: pd.DataFrame,
        bc_refdict: dict,
        p5_refdict: dict,
        p7_refdict: dict | None,
        settings: dict,
    ):

        # Init `rawdict` {idx_seq:{bc_seq:count}}
        self.rawdict = countdict

        # Init `bc_refdict_kseq` {bc_seq:bc_name}
        self.bc_refdict_kseq = {value: key for key, value in bc_refdict.items()}

        # Init `idx_expect_kseq` {idx_seq:idx_name} for expected indexes
        # This will either create p5 only or p5::p7 version
        if p7_refdict is None:
            self.single_index_ref(p5_refdict, runinfo)
        else:
            self.dual_index_ref(p5_refdict, p7_refdict, runinfo)

        # Initialize `idx_known` {idx_name:[]}, empty dict with expected idx_name keys
        self.idx_known = {k: {} for k in self.idx_expect_kseq.values()}

        # Initialize `idx_other` {unknown_idx_seq:{bc_seq:count}}
        self.idx_other = {}

        # Populate `idx_known` and `idx_other` from `rawdict`
        self.combine_known_idx(settings["mismatches"])

        # Generate `bc_set` and `df` from `idx_known`
        self.count_dict_to_dfs()

        # Add run meta to obj and merge runinfo to `df`
        self.add_runinfo(runinfo)

        # Per-index calculations and QC flags
        self.idx_notes = defaultdict(list)
        self.qc_ops(settings)

    #############################
    # init methods
    #############################

    def single_index_ref(self, p5_refdict, runinfo):
        """
        Generate an index reference dictionary for p5-only runs, e/g {'VPX.P5.N':'ACGTACGT'}
        """
        self.idx_expect_kseq = {
            v: k for k, v in p5_refdict.items() if k in runinfo["Barcodes"].values
        }

    def dual_index_ref(self, p5_refdict, p7_refdict, runinfo):
        """
        Generate an index reference dictionary for dual index runs.
        This joins P5::P7 as a single sequence, e/g {'VPX.P5.N_VPX.P7':'ACGTACGTACGTACGT'}
        """
        # Concatenate p7 to p5
        self.idx_expect_kseq = {}

        # Pull only complete rows with these
        cleaned_runinfo = runinfo.dropna(subset=["Barcodes", "(F Barcode)"])

        # Iterate through each row in runinfo to build idx_expect_kseq
        for _, row in cleaned_runinfo.iterrows():
            idx_name = row["Barcodes"]
            p7_name = row["(F Barcode)"]

            # Check if index exists in p5_refdict
            if idx_name not in p5_refdict.keys():
                raise ValueError(f"'{idx_name}' from runinfo not in p5_refdict.")

            idx_seq = p5_refdict[idx_name]  # Get the corresponding index sequence

            # Check if index exists in p7_refdict
            if p7_name not in p7_refdict.keys():
                raise ValueError(f"'{p7_name}' from runinfo not in p7_refdict.")

            p7_seq = p7_refdict[p7_name]  # Get the corresponding p7 sequence

            # Concatenate sequences and create combined name
            full_seq = idx_seq + p7_seq
            full_name = f"{idx_name}_{p7_name}"

            # Add to the dictionary
            self.idx_expect_kseq[full_seq] = full_name

    def combine_known_idx(self, maxdist: int):
        """
        Assign observed index sequences to expected indexes.

        Policy:
        1) Exact match (LD=0) always wins.
        2) Otherwise, assign only if there is a unique best match (smallest LD <= maxdist).
        3) Otherwise (no match or tied best), send to idx_other.
        """
        for obs, counts in self.rawdict.items():
            # 1) exact-match?
            chosen = self.idx_expect_kseq.get(obs)

            if chosen is None:
                # 2) look for unique best match within maxdist
                best_dist = None
                tied = False
                for exp, name in self.idx_expect_kseq.items():
                    d = Levenshtein.distance(obs, exp, score_cutoff=maxdist)
                    if d > maxdist:
                        continue
                    # strictly better match -> adopt it
                    if best_dist is None or d < best_dist:
                        best_dist, chosen, tied = d, name, False
                    # equal best distance -> ambiguous
                    elif d == best_dist:
                        tied = True

                # reject if no match or ambiguous best match
                if best_dist is None or tied:
                    self.idx_other[obs] = counts
                    continue

            # merge counts into chosen index
            dst = self.idx_known.setdefault(chosen, {})
            for bc, n in counts.items():
                dst[bc] = dst.get(bc, 0) + n

        # Sanity check: at least one expected index must have received counts
        if not any(self.idx_known.values()):
            raise ValueError(
                """
                Error! No expected primers found.
                Scroll up to find "primer_path" or "p7_path" for file location.
                Check that runinfo columns match name format in primer reference fasta.
                """
            )

    def count_dict_to_dfs(self):
        """
        Convert dict = {'idx':{'bc':count}} to 3-col df to obtain a per-run set of named barcodes.
            - Uniques are named consistentily across all indexes in rank order.
            - Generates nonredundant df self.bc_set and complete self.named_df with counts.
        """
        # Generate the raw df
        rawdf = pd.DataFrame.from_records(
            [
                (idx_key, bc_key, bc_value)
                for idx_key, idx_dict in self.idx_known.items()
                for bc_key, bc_value in idx_dict.items()
            ],
            columns=["idx_name", "bc_seq", "bc_count"],
        )

        # Dealing with uniques:
        # Get df with total counts per barcode sequence
        df = (
            rawdf.groupby(["bc_seq"])["bc_count"]
            .sum()
            .reset_index()
            .sort_values(by="bc_count", ascending=False)
        )

        # map barcodes names from refdict
        df["bc_name"] = df["bc_seq"].map(self.bc_refdict_kseq)

        # swap unmatched NaN for 'Unique' and add ranks
        df["bc_name"] = df["bc_name"].replace(np.nan, "Unique")
        df = qc.rank_uniques(df=df)

        # Sort by count and select seq/name columns
        df = df.sort_values(by="bc_count", ascending=False)
        df = df[["bc_seq", "bc_name"]]

        # Save a complete set of named barcode seqs for the run
        self.bc_set = df.reset_index(drop=True)

        # Merge df to rawdf to reassign names
        df = pd.merge(rawdf, self.bc_set, how="left", on="bc_seq")
        df = df.sort_values(["idx_name", "bc_count"], ascending=[True, False])

        # If barcode name is found twice, it's in multiple indexes.
        df["x_idx"] = df.duplicated(subset=["bc_name"], keep=False)

        # Set self.named_df
        self.df = df.reset_index(drop=True)

    def add_runinfo(self, runinfo):
        """
        Clean and rename runinfo columns, merge to named_df
        """
        # Get run metadata as object attrs
        self.rnum = runinfo["Run Number"].iloc[0]

        self.rnam = runinfo.loc[
            runinfo.index[runinfo["Run Number"] == "Run Name"].tolist()[0] + 1,
            "Run Number",
        ]

        self.rdat = runinfo.loc[
            runinfo.index[runinfo["Run Number"] == "Date"].tolist()[0] + 1, "Run Number"
        ]

        # Subset and rename per-idx metadata columns from runinfo
        self.idx_meta = runinfo.dropna(subset=["Barcodes", "(F Barcode)"])[
            ["Animal", "Sample", "full_idx", "Date", "Input TOTAL PER BARCODE"]
        ].rename(
            columns={
                "Animal": "samp_group",
                "Sample": "samp_name",
                "full_idx": "idx_name",
                "Date": "samp_date",
                "Input TOTAL PER BARCODE": "input",
            }
        )

        # Merge runinfo to named_df and raise error if any rows are missed
        try:
            merged_df = pd.merge(self.df, self.idx_meta, on="idx_name", how="left")

            # Check if any rows in self.named_df didn't find a match
            if merged_df["samp_group"].isnull().any():
                raise ValueError("rows in named_df did not find idx match in runinfo")

            # Update self.named_df with the merged result
            self.df = merged_df

        except ValueError as e:
            # Re-raise the ValueError with additional context if needed
            raise ValueError(f"Merge operation failed: {str(e)}") from e

    def qc_ops(self, settings: dict):
        """
        Perform QC operations:

        - Global
            - flag x_group barcodes (putative idx_hops)
            - flag short barcodes
            - final column arrangement

        - Index specific
            - calculate proportions
            - flag above-cutoff
            - ldist check for putative parents

        This results in a complete, 'final' long-format dataframe
        experimental: self.idx_notes are built for io.py runinfo updates

        See class & static methods below.
        """
        # flag cross-group barcodes
        self.df["x_group"] = (
            self.df.groupby("bc_name")["samp_group"].transform("nunique") > 1
        )

        # Split df by idx_name to {'idx_name':pd.DataFrame} dict to process individually
        df_dict = {name: group for name, group in self.df.groupby("idx_name")}

        # Per-index processing
        for idx, df in df_dict.items():
            # proportions
            df = self.get_prop(df)

            # above cutoff
            df, note = self.flag_above_cutoff(df)
            note is not None and self.idx_notes[idx].append(note)

            # ldist/putative parents
            df = qc.flag_putative_parents(
                set_name=idx,
                df=df,
                dist=settings["dist_threshold"],
                backend=settings["dist_backend"],
            )

            # collapse ambiguous children
            if settings["collapse"]:
                df, note = self.collapse_ambig_children(df)
                note is not None and self.idx_notes[idx].append(note)

            df_dict[idx] = df

        # Concatenate the modified DataFrames
        self.df = pd.concat(df_dict.values(), axis=0)

        # Flag short barcodes
        self.df = qc.flag_shorts(df=self.df)

        # Arrange columns and rows
        self.df = self.df.loc[
            :,
            [
                "idx_name",  # Index name - P5 or P5_P7
                "samp_group",  # Sample Group - animal, donor, etc
                "samp_name",  # Sample Name - specific sample
                "samp_date",  # From runinfo
                "input",  # ""
                "bc_name",  # Barcode name
                "bc_count",  # Barcode miseq counts
                "proportion",  # Per-index proportion
                "above_cutoff",  # Proportion > 1/input (named) or 100 * min named prop (unique)
                "x_idx",  # Cross-index barcode?
                "x_group",  # Cross-group barcode?
                "short_bc",  # Short barcode?
                "putative_parent",  # ldist hit to named barcode in this index?
                "bc_seq",  # Barcode sequence
            ],
        ].sort_values(["idx_name", "bc_count"], ascending=[True, False])

        # Add run metadata as columns
        add_cols = {"run_num": self.rnum, "run_name": self.rnam, "run_date": self.rdat}

        self.df = pd.concat(
            [pd.DataFrame(add_cols, index=self.df.index), self.df], axis=1
        )

        # Reset index
        self.df.reset_index(drop=True, inplace=True)

    @staticmethod
    def get_prop(df: pd.DataFrame):
        """
        Note: written to run with single-index df via per_index_ops()
        """
        # calculate proportion
        df["proportion"] = df["bc_count"] / df["bc_count"].sum()

        return df

    @staticmethod
    def flag_above_cutoff(df: pd.DataFrame):
        """
        Note: written to run with single-index df via per_index_ops()

        Assign above_cutoff as True for rows where:
            - bc_name does not contain 'Unique_' and proportion is greater than 1/input
            - bc_name contains 'Unique' and bc_count is greater than minimum named proportion * 100
        """
        # Initialize above_cutoff col
        df["above_cutoff"] = False

        # Pull given input
        input_val = df["input"].iloc[0]

        # Collect notes here
        notes = []

        # Handle infinite input gracefully
        if np.isinf(input_val) or np.isnan(input_val):
            named_cutoff = 1 / 10
            notes.append(
                f"Input value is {input_val}, used default cutoff {named_cutoff}"
            )
        else:
            named_cutoff = 1 / input_val

        # Flag above_cutoff for non-unique barcodes
        df.loc[
            ~df["bc_name"].str.contains("Unique.") & (df["proportion"] > named_cutoff),
            "above_cutoff",
        ] = True

        # Get minimum named above_cutoff proportion
        min_named_p = df.loc[df["above_cutoff"], "proportion"].min()

        # Substitute nominal value if min_named_p is missing/nan
        if np.isnan(min_named_p):
            unique_cutoff = named_cutoff
            notes.append(f"No named above_cutoff: used prop > {unique_cutoff:.3g}")
        # Otherwise use 100X minimum named proportion
        else:
            unique_cutoff = 100 * min_named_p

        # Flag above_cutoff for unique barcodes
        df.loc[
            df["bc_name"].str.contains("Unique.") & (df["proportion"] > unique_cutoff),
            "above_cutoff",
        ] = True

        # Join all notes into a single string separated by semicolons
        note = "; ".join(notes) if notes else None

        return df, note

    @staticmethod
    def collapse_ambig_children(df: pd.DataFrame):
        """
        Collapse children with a single 'N' to valid parent rows, and drop multi-'N'.
        Assumes:
        - 'bc_seq' (string with possible 'N's)
        - 'bc_name' (parent name without suffix)
        - 'putative_parent' (now bc_name + '-#' suffix)
        - 'above_cutoff' (bool)
        - 'bc_count' (numeric)
        Returns: (df, note)
        """
        # quick exit if no Ns exist at all
        if not df["bc_seq"].str.contains("N").any():
            return df, None

        note = "ambig_children collapsed"
        df = df.copy()

        # helper: strip a single trailing -<digits> from the putative_parent
        def base_parent_name(s: pd.Series) -> pd.Series:
            # coerce to pandas string dtype (preserves missing as <NA>)
            s2 = s.astype("string")
            # safe replace; .str is valid on string dtype
            return s2.str.replace(r"-.*$", "", regex=True)

        # Drop rows where 'bc_seq' contains more than one 'N'
        df = df[~(df["bc_seq"].str.count("N") > 1)]

        # Identify ambiguous children
        parent_base = base_parent_name(df["putative_parent"].astype("string"))
        eligible_parents = df.loc[df["above_cutoff"], "bc_name"]

        ambig_mask = (
            df["putative_parent"].notna()  # has named putative parent
            & df["bc_seq"].str.contains("N")  # contains a single N (multi was dropped)
            & (~df["above_cutoff"])  # is below cutoff
            & parent_base.isin(eligible_parents)  # parent is above cutoff
        )

        ambig_children = df.loc[ambig_mask]

        # Early exit if nothing to do
        if ambig_children.empty:
            return df, None

        # Remove these plus all remaining rows that contain any N
        df = df[~df["bc_seq"].str.contains("N")]

        # Sum ambiguous children counts by *base* parent name (bc_name without suffix)
        add_counts = (
            ambig_children["bc_count"]
            .groupby(base_parent_name(ambig_children["putative_parent"]))
            .sum()
        )

        # Add those counts to the corresponding parent rows
        df.loc[df["bc_name"].isin(add_counts.index), "bc_count"] += (
            df["bc_name"].map(add_counts).fillna(0)
        )

        # Recompute proportions
        tot = df["bc_count"].sum()
        df["proportion"] = (df["bc_count"] / tot) if tot > 0 else 0.0

        return df, note

    ###################
    # external methods
    ###################

    def samp_group_matrices(self, filt_ac=False):
        """
        Return a per-samp_group dict of barcode proportion matrices,
        optionally filtered for above_cutoff
        """
        df = self.df.copy()

        # Get distinct {samp_group: [idx_name(s)]}
        groups = df.groupby("samp_group")["idx_name"].unique().to_dict()

        # Initialize dict to hold per-samp_group results
        df_mats = {}

        for samp_group, idx_names in groups.items():
            # Optional cutoff filter
            if filt_ac:
                # Pull df copy
                df_filt = df[(df["idx_name"].isin(idx_names))].copy()

                # This is a little overwrought - but I'm leaving it in case we need the flexibility.
                # Conditions for grouping barcodes
                c1 = df_filt["above_cutoff"]  # Above cutoff
                c2 = df_filt["putative_parent"].notna()  # Has putative parent
                c3 = df_filt["bc_name"].str.contains("Unique")  # Is Unique

                # Bitpack the flag columns - this is actually slick.
                df_filt["flag"] = (
                    c1.astype(int) * 4 + c2.astype(int) * 2 + c3.astype(int) * 1
                )

                """
                # Possible flag values:
                flags = {
                    0: 'Below cutoff, No parent, Named',
                    1: 'Below cutoff, No parent, Unique',
                    2: 'Below cutoff, Has parent, Named',
                    3: 'Below cutoff, Has parent, Unique',
                    4: 'Above cutoff, No parent, Named',
                    5: 'Above cutoff, No parent, Unique',
                    6: 'Above cutoff, Has parent, Named',
                    7: 'Above cutoff, Has parent, Unique'
                }

                # With this format, the groups we want are:
                #   0: 4,5 (keep)
                #   1: 0,1,2,3,6,7 (toss)
                """

                # filter above cutoff rows
                df_filt = df_filt[df_filt["flag"].isin([4, 5])].copy()

                # recalculate proportions based on above_cutoff data
                df_filt["proportion"] = df_filt.groupby("idx_name")[
                    "proportion"
                ].transform(lambda x: x / x.sum())

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

            # Get all columns except the first two
            cols_to_sort = df_mat.columns[2:]

            # Sort the DataFrame
            df_mat = df_mat.sort_values(
                by=list(cols_to_sort),
                ascending=[False] * len(cols_to_sort),
                na_position="last",
            ).reset_index(drop=True)

            # Move 'bc_seq' and 'bc_name' to the end
            # df_mat["bc_name"] = df_mat.pop("bc_name")
            # df_mat["bc_seq"] = df_mat.pop("bc_seq")

            # Add to dict
            df_mats[samp_group] = df_mat

        return df_mats

    def samp_group_sidelong_tables(self, legacy_format=False):
        """
        "Sidelong tables" are generated for each samp_group and held in a dict.

        This generates a viewer-oriented table containing adjascent long dataframes
        which indicate cutoffs via padding rows. It's not computer-friendly.

        The legacy format includes a weird header and the original colnames,
        included to ensure continuity with legacy data.
        """
        df = self.df.copy()

        # Get samp_group indexes
        groups = df.groupby("samp_group")["idx_name"].unique().to_dict()

        # Split long df to dict w/ keys = idx_name
        split_dfs = {k: v for k, v in df.groupby("idx_name")}

        # Per-index data formatting
        for idx in split_dfs:
            # subset index's df
            df = split_dfs[idx].copy()

            #########
            # body df
            #########

            df_body = df.copy()[
                [
                    "bc_name",
                    "bc_seq",
                    "bc_count",
                    "proportion",
                    "putative_parent",
                    "x_group",
                    "short_bc",
                    "above_cutoff",
                ]
            ]

            # Conditions for grouping barcodes
            c1 = df_body["above_cutoff"]  # Above cutoff
            c2 = df_body["putative_parent"].notna()  # Has putative parent
            c3 = df_body["bc_name"].str.contains("Unique")  # Is Unique

            # Bitpack the flag columns - this is actually slick.
            df_body["flag"] = (
                c1.astype(int) * 4 + c2.astype(int) * 2 + c3.astype(int) * 1
            )

            """
            # Possible flag values:
            flags = {
                0: 'Below cutoff, No parent, Named',
                1: 'Below cutoff, No parent, Unique',
                2: 'Below cutoff, Has parent, Named',
                3: 'Below cutoff, Has parent, Unique',
                4: 'Above cutoff, No parent, Named',
                5: 'Above cutoff, No parent, Unique',
                6: 'Above cutoff, Has parent, Named',
                7: 'Above cutoff, Has parent, Unique'
            }

            # With this format, the groups we want are:
            #   0: 4,5
            #   1: 6,7
            #   2: 0,1,2,3 (the rest)
            """

            # map the grouping variables
            df_body["rowgroup"] = (
                df_body["flag"].map({4: 0, 5: 0, 6: 1, 7: 1}).fillna(2)
            )

            # sort by ascending group and descending proportion
            df_body = df_body.sort_values(
                by=["rowgroup", "proportion"], ascending=[True, False]
            ).reset_index(drop=True)

            # I'm pained by this:
            # Recalculate group 0's proportions to equal 1...
            ac_tot_prop = df_body.loc[df_body["rowgroup"] == 0, "proportion"].sum()

            if ac_tot_prop != 0:
                df_body.loc[df_body["rowgroup"] == 0, "proportion"] /= ac_tot_prop
            else:
                # Handle edge case if sum is zero
                df_body.loc[df_body["rowgroup"] == 0, "proportion"] = 0

            # find indices where the rowgroup changes
            changes = df_body.index[
                df_body["rowgroup"] != df_body["rowgroup"].shift()
            ].tolist()

            changes.pop(
                0
            )  # artifact of using shift() causes 'change' at first row - pop this.

            # insert empty rows
            pad = pd.DataFrame(
                [[np.nan] * len(df_body.columns)], columns=df_body.columns
            )

            for index in sorted(changes, reverse=True):
                df_body = pd.concat(
                    [df_body.iloc[:index], pad, df_body.iloc[index:]], ignore_index=True
                )

            # drop columns
            df_body = df_body.drop(["flag", "rowgroup", "above_cutoff"], axis=1)

            # Format-specific modification
            if legacy_format:
                # Cast these to string for subsequent map operation
                df_body[["x_group", "short_bc"]] = df_body[
                    ["x_group", "short_bc"]
                ].astype(str)

                # Get mask for spaced rows
                nan_mask = ~df_body["bc_name"].isna()

                # modify (string-cast) boolean columns
                df_body.loc[nan_mask, ["x_group", "short_bc"]] = df_body.loc[
                    nan_mask, ["x_group", "short_bc"]
                ].map(lambda x: "Yes" if x == "1.0" else "" if x == "0.0" else "")

                df_body.loc[~nan_mask, ["x_group", "short_bc"]] = df_body.loc[
                    ~nan_mask, ["x_group", "short_bc"]
                ].map(lambda x: "" if x == "nan" else x)

                # modify putative parent column
                df_body.loc[nan_mask, "putative_parent"] = (
                    df_body.loc[nan_mask, "putative_parent"]
                    .fillna("NA-NA")
                    .map(lambda x: f"{x}" if x != "NA-NA" else x)
                )

                # modify column names
                df_body = df_body.rename(
                    columns={
                        "bc_name": "Barcode",
                        "bc_seq": "Sequence",
                        "bc_count": "Counts",
                        "proportion": "Proportion",
                        "putative_parent": "Hamming Dist to Another Sequence",
                        "x_group": "Multi-group?",
                        "short_bc": "Short barcode?",
                    }
                )

            else:
                # modify boolean columns
                df_body[["x_group", "short_bc"]] = df_body[["x_group", "short_bc"]].map(
                    lambda x: True if x == 1 else False if x == 0 else ""
                )

            # Further padding steps for pd.concat to header:
            # - get column range body + 1
            new_cols = [f"{i}" for i in range(len(df_body.columns) + 1)]

            # - add original colnames as row
            df_body = pd.concat(
                [pd.DataFrame([df_body.columns], columns=df_body.columns), df_body]
            ).reset_index(drop=True)

            # - add padding columns
            for i in range(len(new_cols) - len(df_body.columns)):
                df_body[i] = np.nan

            # - Rename columns for concat
            df_body.columns = new_cols

            ###########
            # header df
            ###########

            # Initialize dict to hold final header values
            header_dict = {}

            # Check these columns for single value
            header_keys = [
                "run_num",  # run-specific
                "run_name",  # ""
                "run_date",  # ""
                "idx_name",  # Index-specific
                "samp_group",  # ""
                "samp_name",  # ""
                "samp_date",  # ""
                "input",  # ""
            ]

            for k in header_keys:
                v = df[k].unique()
                if len(v) != 1:
                    raise ValueError(f"{k} has multiple values: {v}")

                # If unique, update dict
                header_dict[k] = v[0]

            # Format-specific modification
            if legacy_format:
                # Modify date format
                header_dict["samp_date"] = header_dict["samp_date"].strftime("%Y-%m-%d")

                header_data = [
                    # Weird repeated sample metadata
                    [
                        header_dict["samp_group"],
                        header_dict["samp_name"],
                        header_dict["samp_date"],
                        header_dict["samp_group"],
                        header_dict["samp_name"],
                        header_dict["samp_date"],
                    ]
                    + [np.nan] * (len(new_cols) - 6),
                    # Actual "column names"
                    [
                        "Sequencer Input",
                        "Total # Sequences",
                        "# of Barcodes (Excluding Uniques)",
                        "Proportion Uniques",
                    ]
                    + [np.nan] * (len(new_cols) - 4),
                    # Values
                    [
                        header_dict["input"],
                        df["bc_count"].sum(),
                        df[~df["bc_name"].str.contains("Unique")].shape[0],
                        df[df["bc_name"].str.contains("Unique")]["bc_count"].sum()
                        / df["bc_count"].sum(),
                    ]
                    + [np.nan] * (len(new_cols) - 4),
                    # pad
                    [np.nan] * len(new_cols),
                ]

                df_head = pd.DataFrame(header_data, columns=new_cols)

            else:
                # Create the DataFrame with keys and values
                df_head = pd.DataFrame(list(header_dict.items()), columns=new_cols[:2])

                # pad width
                for col in new_cols[2:]:  # Start from the third column of new_cols
                    df_head[col] = np.nan

                # pad length
                pad = pd.DataFrame([[np.nan] * len(new_cols)], columns=new_cols)
                df_head = pd.concat([df_head, pad], ignore_index=True)

                # pad thai?

            # IT'S ALIVE - concat head and body
            df = pd.concat([df_head, df_body], ignore_index=True)

            # Update original dict with final df
            split_dfs[idx] = df

        # Human Centipede - join these sidelong
        sidelongs = {}

        for k, v in groups.items():
            # Extract DataFrames for each key
            dfs = [split_dfs.get(k2) for k2 in v if k2 in split_dfs]

            # Join horizontally and save to new dict
            sidelongs[k] = pd.concat(dfs, axis=1)

        return sidelongs, split_dfs


# %% Versions
"""
v1.6 - moved shared methods w/ xlsx_compiler to utils.qc_ops
v1.61 - Change any "inf" designation to lowercase and pd.to_numeric inputs
v1.7 - Modified parameter passing via settings dict and accomodating changes to flag_putative_parents
v1.8 - Fixed a major bug in combine_known_idx that prevented calls when there
        were close matches in idx_expect (dist 1 between expected caused failed assignment)
v1.9 - Modifed samp_group_matrices to exclude children from above-cutoff matrix
"""
