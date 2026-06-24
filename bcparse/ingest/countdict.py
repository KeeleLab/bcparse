# countdict.py
"""
Name:       countdict.py
Author:     CAG
Version:    1.0
Date:       2026/03/26
Refactored: 2026/05/26

Parse-mode count-dict ingest:
- reference data loading
- count dictionary to `SeqRun` construction
- parse-mode builder state used by workbook emission
"""

# %% Imports

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any, TypeAlias

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein

import bcparse.qc as qc
from bcparse.containers.runinfo import RunInfo
from bcparse.containers.seqrun import SeqRun
from bcparse.containers.seqsamp import SeqSamp

CountDict: TypeAlias = dict[str, dict[str, int]]
RefDict: TypeAlias = dict[str, str]
SettingsDict: TypeAlias = dict[str, Any]

# %% Reference data


def load_parse_reference_data(
    *,
    runinfo_path: str,
    primer_path: str,
    barcode_path: str,
    spike_path: str | None = None,
    p7_path: str | None = None,
) -> tuple[RunInfo, RefDict, RefDict, RefDict | None]:
    print("Formatting reference data... \n")
    runinfo_df = pd.read_excel(runinfo_path)
    runinfo = RunInfo.from_table(
        raw_df=runinfo_df,
        filepath=runinfo_path,
        dual_index=p7_path is not None,
    )
    p5_refdict = _read_csv_to_dict(primer_path)
    bc_refdict = _read_fasta_to_dict(barcode_path)
    if spike_path:
        spike_refdict = _read_fasta_to_dict(spike_path)
        bc_refdict = _merge_reference_dicts(
            primary=bc_refdict,
            extra=spike_refdict,
            extra_label=spike_path,
        )
    p7_refdict = None if p7_path is None else _read_csv_to_dict(p7_path)
    return runinfo, p5_refdict, bc_refdict, p7_refdict


def _read_csv_to_dict(file_path: str) -> dict[str, str]:
    with open(file_path, "r") as f:
        reader = csv.reader(f)
        next(reader)
        return {col1: col2 for col1, col2 in reader}


def _read_fasta_to_dict(file_path: str) -> dict[str, str]:
    sequences, current_header, current_seq = {}, None, []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header, current_seq = line[1:], []
            else:
                current_seq.append(line)
    if current_header is not None:
        sequences[current_header] = "".join(current_seq)
    return sequences


def _merge_reference_dicts(
    *,
    primary: RefDict,
    extra: RefDict,
    extra_label: str,
) -> RefDict:
    overlap = sorted(set(primary).intersection(extra))
    if overlap:
        raise ValueError(
            f"{extra_label!r} contains duplicate barcode name(s): {overlap}"
        )
    return {**primary, **extra}


def build_seqrun_from_countdict(
    *,
    countdict: CountDict,
    runinfo: RunInfo,
    bc_refdict: RefDict,
    p5_refdict: RefDict,
    p7_refdict: RefDict | None,
    settings: SettingsDict,
) -> SeqRun:
    builder = CountDictBuilder(
        countdict=countdict,
        runinfo=runinfo,
        bc_refdict=bc_refdict,
        p5_refdict=p5_refdict,
        p7_refdict=p7_refdict,
        settings=settings,
    ).build()
    if builder.seq_run is None:
        raise RuntimeError("CountDictBuilder did not produce a SeqRun.")
    return builder.seq_run


# %% CountDictBuilder builder class


class CountDictBuilder:
    """
    Internal builder to process the counts dictionary {idx_seq:{bc_seq:count}}
    into a SeqRun plus the legacy intermediate attrs still needed for parse output.
    """

    def __init__(
        self,
        *,
        countdict: CountDict,
        runinfo: RunInfo,
        bc_refdict: RefDict,
        p5_refdict: RefDict,
        p7_refdict: RefDict | None,
        settings: SettingsDict,
    ) -> None:
        # Init `rawdict` {idx_seq:{bc_seq:count}}
        self.rawdict = countdict
        self.runinfo = runinfo
        self.bc_refdict = bc_refdict
        self.p5_refdict = p5_refdict
        self.p7_refdict = p7_refdict
        self.settings = settings
        self.seq_run: SeqRun | None = None
        self.sample_notes: dict[str, list[str]] = defaultdict(list)
        self.sample_meta_df: pd.DataFrame = pd.DataFrame()
        self.bc_refdict_kseq: RefDict = {}
        self.idx_expect_kseq: RefDict = {}
        self.idx_known: CountDict = {}
        self.idx_other: CountDict = {}
        self.bc_set: pd.DataFrame = pd.DataFrame()
        self.df: pd.DataFrame = pd.DataFrame()
        self.rnum: Any = None
        self.rnam: Any = None
        self.rdat: Any = None

    # ====================
    # PARSE
    # ====================

    def build(self) -> "CountDictBuilder":
        """
        Run the full FASTQ/count-dict build path and retain parse-mode state.
        """
        # Init `bc_refdict_kseq` {bc_seq:bc_name}
        self.bc_refdict_kseq = {value: key for key, value in self.bc_refdict.items()}

        # Init `idx_expect_kseq` {idx_seq:idx_name} for expected indexes
        # This will either create p5 only or p5::p7 version
        if self.p7_refdict is None:
            self.single_index_ref(self.p5_refdict, self.runinfo)
        else:
            self.dual_index_ref(self.p5_refdict, self.p7_refdict, self.runinfo)

        # Initialize `idx_known` {idx_name:[]}, empty dict with expected idx_name keys
        self.idx_known = {k: {} for k in self.idx_expect_kseq.values()}

        # Initialize `idx_other` {unknown_idx_seq:{bc_seq:count}}
        self.idx_other = {}

        # Populate `idx_known` and `idx_other` from `rawdict`
        self.combine_known_idx(self.settings["mismatches"])

        # Generate `bc_set` and `df` from `idx_known`
        self.count_dict_to_dfs()

        # Add run meta to obj and merge runinfo to `df`
        self.add_runinfo(self.runinfo)

        # Per-index calculations and QC flags
        self.qc_ops(self.settings)
        self._build_seq_run()

        return self

    def single_index_ref(self, p5_refdict: RefDict, runinfo: RunInfo) -> None:
        """
        Generate an index reference dictionary for p5-only runs, e/g {'VPX.P5.N':'ACGTACGT'}
        """
        self.idx_expect_kseq = {
            v: k
            for k, v in p5_refdict.items()
            if k in runinfo.raw_df["Barcodes"].values
        }

    def dual_index_ref(
        self,
        p5_refdict: RefDict,
        p7_refdict: RefDict,
        runinfo: RunInfo,
    ) -> None:
        """
        Generate an index reference dictionary for dual index runs.
        This joins P5::P7 as a single sequence, e/g {'VPX.P5.N_VPX.P7':'ACGTACGTACGTACGT'}
        """
        # Concatenate p7 to p5
        self.idx_expect_kseq = {}

        # Pull only complete rows with these
        cleaned_runinfo = runinfo.raw_df.dropna(subset=["Barcodes", "(F Barcode)"])

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

    def combine_known_idx(self, maxdist: int) -> None:
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
            if chosen is None:
                raise RuntimeError("Expected index assignment is missing.")
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

    def count_dict_to_dfs(self) -> None:
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

        # Set self.named_df
        self.df = df.reset_index(drop=True)

    def add_runinfo(self, runinfo: RunInfo) -> None:
        """
        Clean normalized runinfo metadata and merge to named_df
        """
        # Get run metadata as object attrs
        self.rnum = runinfo.run_number
        self.rnam = runinfo.run_name
        self.rdat = runinfo.run_date

        # Subset and rename per-idx metadata columns from runinfo
        self.sample_meta_df = runinfo.sample_rows.copy()
        if self.sample_meta_df.empty:
            raise ValueError("RunInfo sample metadata is empty.")

        # Merge runinfo to named_df and raise error if any rows are missed
        try:
            merged_df = pd.merge(
                self.df, self.sample_meta_df, on="idx_name", how="left"
            )

            # Check if any rows in self.named_df didn't find a match
            if merged_df["samp_group"].isnull().any():
                raise ValueError("rows in named_df did not find idx match in runinfo")

            # Update self.named_df with the merged result
            self.df = merged_df
        except ValueError as e:
            raise ValueError(f"Merge operation failed: {str(e)}") from e

    # ====================
    # QC
    # ====================

    def qc_ops(self, settings: SettingsDict) -> None:
        """
        Perform QC operations:

        - Global
            - annotate x_idx as per-sample share of a shared barcode's run total
            - flag short barcodes
            - final column arrangement

        - Index specific
            - calculate proportions
            - flag above-cutoff
            - ldist check for putative parents

        This results in a complete, 'final' long-format dataframe
        experimental: self.sample_notes are built for parse output runinfo updates
        """
        # Split df by idx_name to {'idx_name':pd.DataFrame} dict to process individually
        df_dict = {str(name): group for name, group in self.df.groupby("idx_name")}

        # Per-index processing
        for idx, df in df_dict.items():
            sample_id = df["sample_id"].dropna().iloc[0]
            print(f"{idx} ({sample_id}):", flush=True)

            # proportions
            print(" - calculating proportions", flush=True)
            df = qc.get_prop(df)

            # above cutoff
            print(" - flagging above-cutoff barcodes", flush=True)
            df, note = qc.flag_above_cutoff(df)
            if note is not None:
                self.sample_notes[idx].append(note)

            # ldist/putative parents
            print(" - checking putative parents", flush=True)
            df = qc.flag_putative_parents(
                set_name=idx,
                df=df,
                settings=settings,
                print_header=False,
            )

            # Collapse ambiguous single-N children into eligible parent rows.
            if settings["collapse_ambig"]:
                print(" - collapsing ambiguous children", flush=True)
                df, note = qc.collapse_ambig_children(df)
                if note is not None:
                    self.sample_notes[idx].append(note)

            df_dict[idx] = df

        # Concatenate the modified DataFrames
        self.df = pd.concat(df_dict.values(), axis=0)

        # Flag short barcodes
        self.df = qc.flag_shorts(df=self.df)

        # Add run metadata columns
        add_cols = {
            "run_number": self.rnum,
            "run_name": self.rnam,
            "run_date": self.rdat,
            "filename": Path(settings["sample_path"]).name,
        }
        self.df = pd.concat(
            [pd.DataFrame(add_cols, index=self.df.index), self.df], axis=1
        )

        # For shared barcodes, store each sample's share of that barcode's run total.
        self.df = qc.annotate_x_idx(self.df)

        # Arrange columns and rows
        self.df = self.df.loc[
            :,
            [
                "run_number",
                "run_name",
                "run_date",
                "filename",
                "idx_name",
                "samp_group",
                "samp_name",
                "samp_date",
                "input",
                "bc_name",
                "bc_count",
                "proportion",
                "above_cutoff",
                "x_idx",
                "short_bc",
                "putative_parent",
                "bc_seq",
            ],
        ].sort_values(["idx_name", "bc_count"], ascending=[True, False])

        self.df = qc.normalize_long_df(
            self.df,
            run_number=self.rnum,
            run_name=self.rnam,
            run_date=self.rdat,
            source_path=settings.get("runinfo_path"),
        )

        # Reset index
        self.df.reset_index(drop=True, inplace=True)

    # ====================
    # EMIT
    # ====================

    def _build_seq_run(self) -> None:
        samples: dict[str, SeqSamp] = {}
        for idx_name, sample_df in self.df.groupby("idx_name", sort=False):
            matches = self.runinfo.sample_rows.loc[
                self.runinfo.sample_rows["idx_name"].astype(str).eq(str(idx_name))
            ]
            sample_meta: dict[str, Any]
            if matches.empty:
                sample_meta = {"idx_name": str(idx_name)}
            else:
                sample_meta = {
                    str(key): value for key, value in matches.iloc[0].to_dict().items()
                }
            sample_id = str(sample_df["sample_id"].iloc[0])
            sample_meta["sample_id"] = sample_id
            samples[sample_id] = SeqSamp(
                df=sample_df.reset_index(drop=True),
                sample_meta=sample_meta,
                sample_id=sample_id,
                source="parse",
            )

        self.seq_run = SeqRun(
            samples=samples,
            run_number=self.rnum,
            run_name=self.rnam,
            run_date=self.rdat,
            run_id=self.df["run_id"].iloc[0],
            runinfo=self.runinfo,
            filepath=self.settings["sample_path"],
            source="parse",
            data_df=self.df.copy(),
        )

    def total_reads_by_idx(self) -> dict[str, int]:
        return {
            index: int(sum(barcode.values()))
            for index, barcode in self.idx_known.items()
        }

    def total_barcodes_by_idx(self) -> dict[str, int]:
        return {index: len(barcodes) for index, barcodes in self.idx_known.items()}

    def named_barcodes_by_idx(self) -> dict[str, int]:
        if self.df.empty:
            return {}

        counts = (
            self.df[~self.df["bc_name"].str.contains("Unique", na=False)]
            .groupby("idx_name")
            .size()
        )
        return {str(index): int(count) for index, count in counts.items()}

    def above_cutoff_reads_by_idx(self) -> dict[str, int]:
        if self.df.empty or "above_cutoff" not in self.df.columns:
            return {}

        counts = (
            self.df.loc[self.df["above_cutoff"].fillna(False)]
            .groupby("idx_name")["bc_count"]
            .sum()
        )
        return {str(index): int(count) for index, count in counts.items()}

    def notes_by_idx(self) -> dict[str, str]:
        return {
            index: "; ".join(notes)
            for index, notes in self.sample_notes.items()
            if notes
        }


# %% Versions
"""
v1.0.1 20260623
 - Fixed redundant parse mode logging by adding optional print_header arg flag_putative_parents()

v1.0.0 20260526
 - Initial version, accomodating new model/container structure
 - Supercedes deprecated parse_countdata module and associated code in bcParse.py
 - Refactored into bcparse package structure
"""
