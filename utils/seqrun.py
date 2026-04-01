# seqrun.py

"""
Name:       seqrun.py
Author:     CAG
Version:    1.0
Date:       2026/03/26

`SeqRun` is the main in-memory representation for one sequencing run:
- many `SeqSamp` members
- normalized run-level metadata
- helpers for ingest, parse, and emission

This file includes the core run container plus the parse-mode helper pipeline.
"""

#%% Imports

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from tqdm import tqdm

import utils.qc_ops as qc
from utils.runinfo import RunInfo, make_run_id, make_sample_id
from utils.seqsamp import SeqSamp

#%% SeqRun model

@dataclass
class SeqRun:
    """
    Container for one sequencing run made up of arbitrary samples.
    """
    samples:  dict[str, SeqSamp] = field(default_factory=dict)  # Primary structured representation: one sample per logical sample.
    run_meta: dict[str, Any] = field(default_factory=dict)      # Run-wide labels such as run number, name, date, and filename.
    runinfo:  RunInfo | None = None                             # Optional normalized runinfo metadata companion.
    data_df:  pd.DataFrame | None = None                        # Cached long-format dataframe for callers that already flattened.
    source:   str = "unknown"
    filepath: Path | None = None

    def __post_init__(self) -> None:
        # Make private copies of mutable inputs.
        self.samples = self._coerce_samples(self.samples)
        self.run_meta = dict(self.run_meta)
        self.data_df = None if self.data_df is None else self.data_df.copy()

        if self.runinfo is not None:
            if not isinstance(self.runinfo, RunInfo):
                raise TypeError("SeqRun.runinfo must be a RunInfo instance.")
            if self.filepath is None:
                self.filepath = self.runinfo.filepath
            for key, value in self.runinfo.run_meta.items():
                self.run_meta.setdefault(key, value)

        if self.filepath is not None:
            self.filepath = Path(self.filepath)
            self.run_meta.setdefault("filepath", self.filepath)
            self.run_meta.setdefault("filename", self.filepath.name)

        self.run_meta.setdefault(
            "run_id",
            make_run_id(
                self.run_meta.get("run_number"),
                self.run_meta.get("filename"),
                self.run_meta.get("run_name"),
            ),
        )

    # ====================
    # INGEST
    # ====================

    @classmethod
    def from_long_df(
        cls,
        df: pd.DataFrame,
        *,
        source: str = "unknown",
        runinfo: RunInfo | None = None,
    ) -> "SeqRun":
        """
        Materialize a SeqRun from a normalized long-format dataframe representing one run.
        """
        df = qc.normalize_long_df(df)

        if df.empty:
            return cls(source=source)

        group_col = cls.long_df_run_group_col(df)
        if group_col is not None:
            run_ids = df[group_col].dropna().astype(str).unique().tolist()
            if len(run_ids) > 1:
                raise ValueError(
                    f"SeqRun.from_long_df expected one {group_col}, found {run_ids}"
                )

        run_meta = {}
        for col in ["run_id", "run_number", "run_name", "run_date", "filename"]:
            if col in df.columns:
                # Long-format metadata repeats across rows, so take the first real value.
                non_null = df[col].dropna()
                value = None if non_null.empty else non_null.iloc[0]
                if value is not None:
                    run_meta[col] = value

        for col in ["run_id", "run_name", "run_date", "filename"]:
            if col not in df.columns:
                continue
            non_null = df[col].dropna().astype(str).unique().tolist()
            if len(non_null) > 1:
                raise ValueError(
                    f"SeqRun.from_long_df expected one {col}, found {non_null}"
                )

        samples = {}
        for sample_id, sample_df in df.groupby("sample_id", dropna=False, sort=False):
            sample_meta = {}
            for col in [
                "sample_id",
                "idx_name",
                "samp_group",
                "samp_name",
                "samp_date",
                "input",
                "run_number",
                "run_name",
                "run_date",
                "filename",
            ]:
                if col not in sample_df.columns:
                    continue
                # Sample metadata repeats across rows, so take the first real value.
                non_null = sample_df[col].dropna()
                value = None if non_null.empty else non_null.iloc[0]
                if value is not None:
                    sample_meta[col] = value

            resolved_sample_id = None if pd.isna(sample_id) else str(sample_id)
            samples[resolved_sample_id or sample_meta.get("sample_id") or str(len(samples))] = SeqSamp.from_table(
                sample_df=sample_df,
                sample_meta=sample_meta,
                sample_id=resolved_sample_id,
                source=source,
            )

        return cls(
            samples=samples,
            run_meta=run_meta,
            runinfo=runinfo,
            source=source,
        )

    @staticmethod
    def _coerce_samples(
        samples: Mapping[str, SeqSamp] | Iterable[SeqSamp],
    ) -> dict[str, SeqSamp]:
        """
        Normalize incoming sample containers to the internal `{sample_id: SeqSamp}` shape.
        """
        if isinstance(samples, Mapping):
            items = samples.items()
        else:
            items = ((None, sample) for sample in samples)

        sample_map: dict[str, SeqSamp] = {}
        for key, sample in items:
            resolved_id = key or sample.sample_id or sample.sample_meta.get("sample_id")
            if resolved_id is None:
                raise ValueError("SeqRun samples require a resolvable sample_id.")
            resolved_id = str(resolved_id)
            if resolved_id in sample_map:
                raise ValueError(f"Duplicate sample_id in SeqRun: {resolved_id}")

            if sample.sample_id != resolved_id:
                sample = SeqSamp(
                    df=sample.df,
                    sample_meta=sample.sample_meta,
                    sample_id=resolved_id,
                    source=sample.source,
                )

            sample_map[resolved_id] = sample

        return sample_map

    @classmethod
    def long_df_run_group_col(cls, df: pd.DataFrame) -> str | None:
        """
        Pick the best available column for splitting long-format data into runs.
        """
        df = qc.normalize_long_df(df)

        for col in ["run_id", "run_number", "filename", "run_name"]:
            if col in df.columns and df[col].notna().any():
                return col

        return None

    @staticmethod
    def annotate_x_idx(df: pd.DataFrame) -> pd.DataFrame:
        """
        For shared barcodes, store each sample's share of that barcode's run total.

        Barcodes confined to one sample are left blank (`NaN`).
        """
        df = qc.normalize_long_df(df)
        if df.empty:
            return df

        df = df.copy()
        group_cols = ["run_number", "bc_name"]
        n_samples = df.groupby(group_cols)["idx_name"].transform("nunique")
        bc_total = df.groupby(group_cols)["bc_count"].transform("sum")

        df["x_idx"] = np.where(
            n_samples > 1,
            df["bc_count"] / bc_total.replace(0, np.nan),
            np.nan,
        )

        return df

    # ====================
    # PARSE
    # ====================

    @classmethod
    def load_parse_reference_data(
        cls,
        *,
        runinfo_path: str,
        primer_path: str,
        barcode_path: str,
        spike_path: str | None = None,
        p7_path: str | None = None,
    ) -> tuple[RunInfo, dict, dict, dict | None]:
        """
        Load parse-mode reference metadata, primers, and barcode names from disk.
        """
        print("Formatting reference data... \n")

        runinfo_df = pd.read_excel(runinfo_path)
        runinfo = RunInfo.from_table(
            raw_runinfo=runinfo_df,
            filepath=runinfo_path,
            dual_index=p7_path is not None,
        )
        p5_refdict = cls._read_csv_to_dict(primer_path)
        bc_refdict = cls._read_fasta_to_dict(barcode_path)
        if spike_path:
            spike_refdict = cls._read_fasta_to_dict(spike_path)
            bc_refdict = cls._merge_reference_dicts(
                primary=bc_refdict,
                extra=spike_refdict,
                extra_label=spike_path,
            )
        p7_refdict = None if p7_path is None else cls._read_csv_to_dict(p7_path)

        return runinfo, p5_refdict, bc_refdict, p7_refdict

    @classmethod
    def parse_fastq(
        cls,
        *,
        countdict: dict,
        runinfo: RunInfo,
        bc_refdict: dict,
        p5_refdict: dict,
        p7_refdict: dict | None,
        settings: dict,
    ) -> "_ParseRunFastq":
        """
        Build the parse-mode helper that owns FASTQ-specific intermediate state.
        """
        return _ParseRunFastq(
            countdict=countdict,
            runinfo=runinfo,
            bc_refdict=bc_refdict,
            p5_refdict=p5_refdict,
            p7_refdict=p7_refdict,
            settings=settings,
        ).build()

    @staticmethod
    def _read_csv_to_dict(file_path: str) -> dict[str, str]:
        """
        Parse a two-column csv into a sequence/name dict.
        """
        result = {}
        with open(file_path, "r") as csvfile:
            reader = csv.reader(csvfile)
            next(reader)
            for row in reader:
                col1, col2 = row
                result[col1] = col2
        return result

    @staticmethod
    def _read_fasta_to_dict(file_path: str) -> dict[str, str]:
        """
        Parse a standard fasta into a header/sequence dict.
        """
        sequences = {}
        current_seq = []
        current_header = None

        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences[current_header] = "".join(current_seq)
                        current_seq = []
                    current_header = line[1:]
                else:
                    current_seq.append(line)

            if current_header is not None:
                sequences[current_header] = "".join(current_seq)

        return sequences

    @staticmethod
    def _merge_reference_dicts(
        *,
        primary: dict[str, str],
        extra: dict[str, str],
        extra_label: str,
    ) -> dict[str, str]:
        """
        Merge an additional fasta-derived reference dict into the primary set.
        """
        overlap = sorted(set(primary).intersection(extra))
        if overlap:
            raise ValueError(
                f"Additional reference {extra_label!r} contains duplicate barcode name(s): {overlap}"
            )

        merged = primary.copy()
        merged.update(extra)
        return merged

    # ====================
    # EMIT
    # ====================

    @property
    def sample_list(self) -> list[SeqSamp]:
        return list(self.samples.values())

    @property
    def sample_ids(self) -> list[str]:
        return list(self.samples.keys())

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    @property
    def sample_meta_df(self) -> pd.DataFrame:
        if not self.samples:
            return pd.DataFrame()
        return pd.DataFrame([sample.sample_meta for sample in self.sample_list])

    @property
    def df(self) -> pd.DataFrame:
        """
        Flatten the run back to one long-format dataframe.
        """
        if self.data_df is not None:
            return self.data_df.copy()
        if not self.samples:
            return pd.DataFrame()
        df = pd.concat([sample.df for sample in self.sample_list], ignore_index=True)
        for key, value in self.run_meta.items():
            if key not in df.columns or df[key].isna().all():
                df[key] = value
        return df

    @property
    def filename(self) -> str | None:
        value = self.run_meta.get("filename")
        return None if value is None else str(value)

    @property
    def run_number(self) -> Any:
        return self.run_meta.get("run_number")

    @property
    def run_name(self) -> Any:
        return self.run_meta.get("run_name")

    @property
    def run_date(self) -> Any:
        return self.run_meta.get("run_date")

    @property
    def df_above_cutoff(self) -> pd.DataFrame:
        """
        Convenience view of only above-cutoff rows.
        """
        df = self.df
        if "above_cutoff" not in df.columns:
            return df.iloc[0:0].copy()
        return df.loc[df["above_cutoff"].fillna(False)].copy()

    @property
    def groups(self) -> list[str]:
        """
        Return distinct sample groups represented in this run.
        """
        df = self.df
        if "samp_group" not in df.columns:
            return []
        return df["samp_group"].dropna().astype(str).unique().tolist()

    def get_metadata(self, key: str, default: Any = None) -> Any:
        return self.run_meta.get(key, default)

    def copy(self, **updates: Any) -> "SeqRun":
        """
        Deep-ish copy used when callers want another mutable run object.
        """
        data = {
            "samples": {
                sample_id: SeqSamp(
                    df=sample.df.copy(),
                    sample_meta=sample.sample_meta.copy(),
                    sample_id=sample.sample_id,
                    source=sample.source,
                )
                for sample_id, sample in self.samples.items()
            },
            "run_meta": self.run_meta.copy(),
            "runinfo": None
            if self.runinfo is None
            else RunInfo(
                run_meta=self.runinfo.run_meta.copy(),
                sample_meta={k: v.copy() for k, v in self.runinfo.sample_meta.items()},
                raw_df=None if self.runinfo.raw_df is None else self.runinfo.raw_df.copy(),
                source=self.runinfo.source,
                filepath=self.runinfo.filepath,
            ),
            "data_df": None if self.data_df is None else self.data_df.copy(),
            "source": self.source,
            "filepath": self.filepath,
        }
        data.update(updates)
        return SeqRun(**data)

#%% _ParseRunFastq builder class

class _ParseRunFastq:
    """
    Internal builder to process the counts dictionary {idx_seq:{bc_seq:count}}
    into a SeqRun plus the legacy intermediate attrs still needed for parse output.
    """

    def __init__(
        self,
        *,
        countdict: dict,
        runinfo: RunInfo,
        bc_refdict: dict,
        p5_refdict: dict,
        p7_refdict: dict | None,
        settings: dict,
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
        self.sample_meta_df: pd.DataFrame | None = None

    # ====================
    # PARSE
    # ====================

    def build(self) -> "_ParseRunFastq":
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

    def single_index_ref(self, p5_refdict, runinfo):
        """
        Generate an index reference dictionary for p5-only runs, e/g {'VPX.P5.N':'ACGTACGT'}
        """
        self.idx_expect_kseq = {
            v: k for k, v in p5_refdict.items() if k in runinfo.raw_df["Barcodes"].values
        }

    def dual_index_ref(self, p5_refdict, p7_refdict, runinfo):
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

        # Set self.named_df
        self.df = df.reset_index(drop=True)

    def add_runinfo(self, runinfo):
        """
        Clean normalized runinfo metadata and merge to named_df
        """
        # Get run metadata as object attrs
        self.rnum = runinfo.run_meta.get("run_number")
        self.rnam = runinfo.run_meta.get("run_name")
        self.rdat = runinfo.run_meta.get("run_date")

        # Subset and rename per-idx metadata columns from runinfo
        self.sample_meta_df = runinfo.sample_meta_df.copy()
        if self.sample_meta_df.empty:
            raise ValueError("RunInfo sample metadata is empty.")

        # Merge runinfo to named_df and raise error if any rows are missed
        try:
            merged_df = pd.merge(self.df, self.sample_meta_df, on="idx_name", how="left")

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

    def qc_ops(self, settings: dict):
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
        df_dict = {name: group for name, group in self.df.groupby("idx_name")}

        # Per-index processing
        for idx, df in df_dict.items():
            with tqdm(total=1, bar_format="{desc} | {elapsed}", leave=True) as pbar:
                pbar.set_description(f"Processing {idx}")

            # proportions
            df = qc.get_prop(df)

            # above cutoff
            df, note = qc.flag_above_cutoff(df)
            note is not None and self.sample_notes[idx].append(note)

            # ldist/putative parents
            df = qc.flag_putative_parents(
                set_name=idx,
                df=df,
                settings=settings,
            )

            # collapse ambiguous children
            if settings["collapse"]:
                df, note = qc.collapse_ambig_children(df)
                note is not None and self.sample_notes[idx].append(note)

            df_dict[idx] = df

        # Concatenate the modified DataFrames
        self.df = pd.concat(df_dict.values(), axis=0)

        # Flag short barcodes
        self.df = qc.flag_shorts(df=self.df)

        # Add run metadata columns
        add_cols = {"run_number": self.rnum, "run_name": self.rnam, "run_date": self.rdat}
        self.df = pd.concat([pd.DataFrame(add_cols, index=self.df.index), self.df], axis=1)

        # For shared barcodes, store each sample's share of that barcode's run total.
        self.df = SeqRun.annotate_x_idx(self.df)

        # Arrange columns and rows
        self.df = self.df.loc[
            :,
            [
                "run_number",
                "run_name",
                "run_date",
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
        )

        # Parse-mode identity is idx_name, not the compile-mode composite sample_id.
        self.df["sample_id"] = self.df["idx_name"]

        # Reset index
        self.df.reset_index(drop=True, inplace=True)

    # ====================
    # EMIT
    # ====================

    def _build_seq_run(self):
        run_meta = {
            "run_number": self.rnum,
            "run_name": self.rnam,
            "run_date": self.rdat,
        }

        samples = {}
        for idx_name, sample_df in self.df.groupby("idx_name", sort=False):
            sample_meta = dict(
                self.runinfo.get_sample_meta_by_idx_name(idx_name, {"idx_name": idx_name})
            )
            sample_id = idx_name
            sample_meta["sample_id"] = sample_id
            samples[sample_id] = SeqSamp(
                df=sample_df.reset_index(drop=True),
                sample_meta=sample_meta,
                sample_id=sample_id,
                source="parse",
            )

        self.seq_run = SeqRun(
            samples=samples,
            run_meta=run_meta,
            runinfo=self.runinfo,
            source="parse",
            data_df=self.df.copy(),
        )

    def total_reads_by_idx(self) -> dict[str, int]:
        return {
            index: sum(barcode.values())
            for index, barcode in self.idx_known.items()
        }

    def total_barcodes_by_idx(self) -> dict[str, int]:
        return {
            index: len(barcodes)
            for index, barcodes in self.idx_known.items()
        }

    def named_barcodes_by_idx(self) -> dict[str, int]:
        if self.df.empty:
            return {}

        return (
            self.df[~self.df["bc_name"].str.contains("Unique", na=False)]
            .groupby("idx_name")
            .size()
            .to_dict()
        )

    def above_cutoff_reads_by_idx(self) -> dict[str, int]:
        if self.df.empty or "above_cutoff" not in self.df.columns:
            return {}

        return (
            self.df.loc[self.df["above_cutoff"].fillna(False)]
            .groupby("idx_name")["bc_count"]
            .sum()
            .to_dict()
        )

    def notes_by_idx(self) -> dict[str, str]:
        return {
            index: "; ".join(notes)
            for index, notes in self.sample_notes.items()
            if notes
        }


    def write_output(
        self,
        *,
        settings: dict,
        runinfo_df: pd.DataFrame | None = None,
    ) -> None:
        """
        Write parse-mode csv/xlsx outputs from this parsed FASTQ run.
        """
        from utils.data_views import samp_group_matrices, sidelong_tables

        out_path = settings["out_path"]
        legacy_format = settings["legacy_format"]
        filt_mat_ac = settings["filt_mat_ac"]
        collapse_to_parent = settings["collapse_to_parent"]
        seq_run = self.seq_run
        runinfo_obj = self.runinfo

        if seq_run is None:
            raise ValueError("Parse output requested before SeqRun was built.")

        # Initialize runinfo output table.
        if runinfo_df is None:
            runinfo_df = pd.DataFrame() if runinfo_obj is None else runinfo_obj.copy_raw_df()
        runinfo_df = qc.format_output_dates(runinfo_df, columns=["Date", "run_date", "samp_date"])

        os.makedirs(out_path, exist_ok=True)

        # Write full long df for downstream re-use / inspection.
        qc.format_output_dates(
            seq_run.df,
            columns=["run_date", "samp_date"],
        ).to_csv(os.path.join(out_path, f"{seq_run.run_name}_concat.csv"), index=False)

        # Build matrix and sample-table workbook views.
        group_mats_dict = samp_group_matrices(
            seq_run,
            filt_ac=filt_mat_ac,
            collapse_to_parent=collapse_to_parent,
        )

        if collapse_to_parent:
            group_tabs_dict, _ = sidelong_tables(
                pd.concat(
                    [sample.collapse_to_parent().df for sample in seq_run.sample_list],
                    ignore_index=True,
                ),
                legacy_format=legacy_format,
            )
            group_tabs_raw_dict, _ = sidelong_tables(
                seq_run.df.copy(),
                legacy_format=legacy_format,
            )
        else:
            group_tabs_dict, _ = sidelong_tables(
                seq_run.df.copy(),
                legacy_format=legacy_format,
            )
            group_tabs_raw_dict = None

        expected_groups = [] if runinfo_obj is None else runinfo_obj.expected_groups
        groups = seq_run.groups
        missing = set(expected_groups) - set(groups)

        if missing:
            print(
                f""" 
        No barcodes found for {missing}!
        """
            )

        if not set(group_mats_dict.keys()) == set(groups) or not set(group_tabs_dict.keys()) == set(groups):
            raise ValueError(
                f"""
        Check concat output for missing data. 
        Details:
        - group_mats_dict keys: {set(group_mats_dict.keys())}
        - group_tabs_dict keys: {set(group_tabs_dict.keys())}
        - expected groups: {set(groups)}

        Previously seen issues here:
        - ran 239M/M2 instead of 239M/M2_dual-index.
        """
            )

        # Add analysis-specific summary columns to runinfo.
        runinfo_df = self._modify_runinfo(runinfo_df, seq_run)

        with pd.ExcelWriter(os.path.join(out_path, f"{seq_run.run_name}_Analysis.xlsx")) as writer:
            runinfo_df.to_excel(writer, sheet_name="runinfo", index=False)

            for samp_group in groups:
                print(f"Writing {samp_group}...")

                if group_tabs_raw_dict is not None:
                    group_tabs_raw_dict[samp_group].to_excel(
                        writer, sheet_name=f"{samp_group} raw", index=False, header=False
                    )

                group_tabs_dict[samp_group].to_excel(
                    writer, sheet_name=f"{samp_group} samples", index=False, header=False
                )

                group_mats_dict[samp_group].to_excel(
                    writer, sheet_name=f"{samp_group} matrix", index=True, header=True
                )

            settings["analysis date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            settings_df = pd.DataFrame(list(settings.items()))
            settings_df.to_excel(writer, sheet_name="Analysis settings", index=False, header=False)


    def _modify_runinfo(
        self,
        runinfo_df: pd.DataFrame,
        seq_run: SeqRun,
    ) -> pd.DataFrame:
        """
        Add parse summary columns and notes to the output runinfo table.
        """
        modified_runinfo = runinfo_df.copy()
        runinfo_obj = seq_run.runinfo
        full_index_col = "full_idx" if runinfo_obj is None else runinfo_obj.full_index_col
        idx_notes = self.notes_by_idx()

        modified_runinfo["n_reads"] = modified_runinfo[full_index_col].map(
            self.total_reads_by_idx()
        )
        modified_runinfo["n_bc"] = modified_runinfo[full_index_col].map(
            self.total_barcodes_by_idx()
        )
        modified_runinfo["n_bc_named"] = modified_runinfo[full_index_col].map(
            self.named_barcodes_by_idx()
        )
        modified_runinfo["n_reads_ac"] = modified_runinfo[full_index_col].map(
            self.above_cutoff_reads_by_idx()
        )
        modified_runinfo["count-to-input_ratio"] = (
            modified_runinfo["n_reads"] / modified_runinfo["Input TOTAL PER BARCODE"]
        )

        if idx_notes:
            modified_runinfo["notes"] = modified_runinfo[full_index_col].map(idx_notes)

        return modified_runinfo

#%% Versions
"""
v1.0 
- Initial version, accomodating new model/container structure
- Supercedes deprecated parse_countdata module and associated code in bcParse.py
"""
