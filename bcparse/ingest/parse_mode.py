# parse_mode.py
"""
Name:       parse_mode.py
Author:     CAG
Version:    1.0.0
Date:       2026/06/25

Parse-mode pipeline for one FASTQ plus one runinfo workbook.

This module is intentionally linear, ParseMode.__init__ describes action.

The goal is to keep parse-mode behavior discoverable in one place without
turning the light container classes into alternate pipeline entry points.

Major classes and helpers:

    ParseMode
        Owns settings expansion, reference loading, FASTQ counting, QC, and
        final SeqRun construction for normal parse-mode execution.

    stream_fastq / stream_to_countdict
        Read FASTQ records and collect observed index/barcode counts.

    parse_fastq_rec
        Extracts the target index and barcode sequences from one FASTQ record.
"""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from importlib.metadata import version as package_version
from pathlib import Path
from typing import Any, Generator, Tuple, TypeAlias, cast

import numpy as np
import pandas as pd
import regex
from rapidfuzz.distance import Levenshtein
from tqdm import tqdm

import bcparse.qc as qc
from bcparse.config import get_stock_settings
from bcparse.containers.runinfo import RunInfo
from bcparse.containers.seqrun import SeqRun
from bcparse.containers.seqsamp import SeqSamp
from bcparse.settings import ParseSettings

CountDict: TypeAlias = dict[str, dict[str, int]]
RefDict: TypeAlias = dict[str, str]
SettingsDict: TypeAlias = dict[str, Any]

#%% ParseMode

class ParseMode:
    """
    Top-level parse-mode orchestration for one sequencing run.

    Call site:
        bcparse.__main__._run_parse() builds this class as:
            parse_run = ParseMode(settings=parse_settings)

    Inputs:
        ParseSettings object from CLI/user configuration:
            - sample FASTQ path
            - runinfo workbook path
            - stock name
            - dual-index
            - flag, spike, and parse/QC thresholds.

    Outputs:
        seq_run:
            Final SeqRun container with SeqSamp children, RunInfo metadata, the
            normalized long dataframe, source='parse', and filepath set to the
            parsed FASTQ.
        df:
            Final normalized long-format dataframe used by SeqRun and workbook
            emission.
        summary helpers:
            total_reads_by_idx(), total_barcodes_by_idx(),
            named_barcodes_by_idx(), above_cutoff_reads_by_idx(), and
            notes_by_idx() provide runinfo/workbook summary columns.
    """

    settings: ParseSettings  # Original user/CLI parse config.
    runtime_settings: SettingsDict  # Stock-expanded runtime config.
    sample_notes: dict[str, list[str]]  # QC notes keyed by idx_name.

    runinfo: RunInfo  # Validated run/sample metadata.
    p5_refdict: RefDict  # P5 name -> sequence reference.
    bc_refdict: RefDict  # Barcode name -> sequence reference.
    p7_refdict: RefDict | None  # Optional P7 name -> sequence reference.
    bc_refdict_kseq: RefDict  # Barcode sequence -> name lookup.

    raw_countdict: CountDict  # Observed idx_seq -> bc_seq counts.
    idx_expect_kseq: RefDict  # Expected idx_seq -> idx_name lookup.
    idx_known: CountDict  # Assigned idx_name -> bc_seq counts.
    idx_other: CountDict  # Unassigned idx_seq -> bc_seq counts.

    bc_set: pd.DataFrame  # Nonredundant run barcode seq/name table.
    sample_meta_df: pd.DataFrame  # Normalized RunInfo sample metadata.
    df: pd.DataFrame  # Final normalized long-format parse dataframe.

    rnum: Any  # Run number copied from RunInfo.
    rnam: Any  # Run name copied from RunInfo.
    rdat: Any  # Run date copied from RunInfo.

    seq_run: SeqRun  # Final parse-mode SeqRun container.

    def __init__(self, *, settings: ParseSettings) -> None:

        self.settings = settings
        self.sample_notes = defaultdict(list)

        # 1. Expand ParseSettings into runtime parse configuration.
        self.runtime_settings = self._build_runtime_settings(settings)

        # 2. Load runinfo and reference data.
        self.p5_refdict, self.bc_refdict, self.p7_refdict = self._load_reference_data()

        # 3. Validate runinfo before counting.
        self.runinfo = self._load_runinfo()
        self._validate_runinfo()

        # 4. Stream FASTQ into a raw count dictionary.
        print("\n")
        self.raw_countdict = stream_to_countdict(settings=self.runtime_settings)
        print("\n")

        # 5. Build index references and assign counts.
        self._build_index_references()
        self._assign_observed_indexes(maxdist=int(self.runtime_settings["mismatches"]))

        # 6. Convert count dictionary data into a long barcode/count dataframe.
        self._count_dict_to_dfs()

        # 8. Merge normalized runinfo sample metadata.
        self._add_runinfo()

        # 9. Run per-sample QC.
        self.qc_ops(self.runtime_settings)

        # 10. Build seq_run container.
        self._build_seq_run()


    def _build_runtime_settings(self, settings: ParseSettings) -> dict:

        settings.validate()

        runtime_settings = get_stock_settings(
            settings.stock,
            settings.dualindex,
            append_spike_ref=settings.append_spike_ref,
        )

        runtime_settings.update(settings.to_dict())

        runtime_settings["software_version"] = package_version("bcparse")
        return runtime_settings


    def _load_reference_data(self) -> tuple[RefDict, RefDict, RefDict | None]:

        print("Formatting reference data... \n")

        p7_path = (
            str(self.runtime_settings["p7_path"]) if self.settings.dualindex else None
        )

        p5_refdict = self._read_csv_to_dict(str(self.runtime_settings["primer_path"]))
        bc_refdict = self._read_fasta_to_dict(
            str(self.runtime_settings["barcode_path"])
        )
        spike_path = (
            self.runtime_settings.get("spike_path")
            if self.settings.append_spike_ref
            else None
        )
        if spike_path:
            spike_refdict = self._read_fasta_to_dict(str(spike_path))
            bc_refdict = self._merge_reference_dicts(
                primary=bc_refdict,
                extra=spike_refdict,
                extra_label=str(spike_path),
            )

        p7_refdict = None if p7_path is None else self._read_csv_to_dict(p7_path)
        return p5_refdict, bc_refdict, p7_refdict


    def _read_csv_to_dict(self, file_path: str) -> RefDict:
        with open(file_path, "r") as f:
            reader = csv.reader(f)
            next(reader)
            return {col1: col2 for col1, col2 in reader}


    def _read_fasta_to_dict(self, file_path: str) -> RefDict:
        sequences: RefDict = {}
        current_header: str | None = None
        current_seq: list[str] = []
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
        self,
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


    def _load_runinfo(self) -> RunInfo:
        runinfo_df = pd.read_excel(self.settings.runinfo_path)
        return RunInfo.from_table(
            raw_df=runinfo_df,
            filepath=self.settings.runinfo_path,
            dual_index=self.settings.dualindex,
        )


    def _validate_runinfo(self) -> None:
        collisions = self.runinfo.find_idx_sample_id_collisions()
        if not collisions:
            return

        lines = [
            "Error: duplicate sample identities found in runinfo.",
            "Distinct idx_name values resolve to the same sample_id:",
            "",
        ]
        for sample_id, idx_names in collisions.items():
            lines.append(f"- {sample_id}")
            lines.append(f"  idx_name(s): {', '.join(idx_names)}")

        raise ValueError("\n".join(lines))


    def _build_index_references(self) -> None:

        # Init `bc_refdict_kseq` {bc_seq:bc_name}
        self.bc_refdict_kseq = {value: key for key, value in self.bc_refdict.items()}

        # Init `idx_expect_kseq` {idx_seq:idx_name} for expected indexes
        # This will either create p5 only or p5::p7 version
        self.idx_expect_kseq = {}
        if self.p7_refdict is None:
            self._single_index_ref()
        else:
            self._dual_index_ref()

        # Init `idx_known` {idx_name:[]}, empty dict with expected idx_name keys
        self.idx_known = {k: {} for k in self.idx_expect_kseq.values()}

        # Initialize `idx_other` {unknown_idx_seq:{bc_seq:count}}
        self.idx_other = {}


    def _single_index_ref(self) -> None:
        """
        Generate an index reference dictionary for p5-only runs, e/g {'VPX.P5.N':'ACGTACGT'}
        """
        self.idx_expect_kseq = {
            v: k
            for k, v in self.p5_refdict.items()
            if k in self.runinfo.raw_df["Barcodes"].values
        }


    def _dual_index_ref(self) -> None:
        """
        Generate an index reference dictionary for dual index runs.
        This joins P5::P7 as a single sequence, e/g {'VPX.P5.N_VPX.P7':'ACGTACGTACGTACGT'}
        """
        p7_refdict = self.p7_refdict
        if p7_refdict is None:
            raise ValueError("Dual-index parse requires a P7 reference.")

        # Concatenate p7 to p5
        self.idx_expect_kseq = {}

        # Pull only complete rows with these
        cleaned_runinfo = self.runinfo.raw_df.dropna(subset=["Barcodes", "(F Barcode)"])

        # Iterate through each row in runinfo to build idx_expect_kseq
        for _, row in cleaned_runinfo.iterrows():
            idx_name = row["Barcodes"]
            p7_name = row["(F Barcode)"]

            # Check if index exists in p5_refdict
            if idx_name not in self.p5_refdict.keys():
                raise ValueError(f"'{idx_name}' from runinfo not in p5_refdict.")

            idx_seq = self.p5_refdict[idx_name]  # Get the corresponding index sequence

            # Check if index exists in p7_refdict
            if p7_name not in p7_refdict.keys():
                raise ValueError(f"'{p7_name}' from runinfo not in p7_refdict.")

            p7_seq = p7_refdict[p7_name]  # Get the corresponding p7 sequence

            # Concatenate sequences and create combined name
            full_seq = idx_seq + p7_seq
            full_name = f"{idx_name}_{p7_name}"

            # Add to the dictionary
            self.idx_expect_kseq[full_seq] = full_name


    def _assign_observed_indexes(self, maxdist: int) -> None:
        """
        Assign observed index sequences to expected indexes.

        Policy:
        1) Exact match (LD=0) always wins.
        2) Otherwise, assign only if there is a unique best match (smallest LD <= maxdist).
        3) Otherwise (no match or tied best), send to idx_other.
        """
        for obs, counts in self.raw_countdict.items():
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


    def _count_dict_to_dfs(self) -> None:
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


    def _add_runinfo(self) -> None:

        # Get run metadata as object attrs
        self.rnum = self.runinfo.run_number
        self.rnam = self.runinfo.run_name
        self.rdat = self.runinfo.run_date

        # Subset and rename per-idx metadata columns from runinfo
        self.sample_meta_df = self.runinfo.sample_rows.copy()
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
            filepath=self.settings.sample_path,
            source="parse",
            data_df=self.df.copy(),
        )

    # Non-init callable summaries, used in emit

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

# %% Fastq handling

def stream_to_countdict(settings: dict):
    """
    Use stream_fastq() and parse_fastq_rec class to collect {idx:{bc:count}}
    """
    res_pfq = {}

    # Initialize search patterns
    parse_fastq_rec.init_patterns(settings)

    for _, seq, _, qual in stream_fastq(settings["sample_path"]):
        res = parse_fastq_rec(
            seq=seq,
            qual=qual,
            settings=settings,
        )

        # Ignore read if barcode is missing or low mean q
        if not hasattr(res, "bc_seq") or res.bc_qual < settings["mean_qual"]:
            continue

        # Process p5 - skip if missing or low mean q
        if not (hasattr(res, "p5_seq") and hasattr(res, "p5_qual")):
            continue

        p5_seq = res.p5_seq if res.p5_qual >= settings["mean_qual"] else None
        if p5_seq is None:
            continue

        # Optionally process p7, fill instead of skip:
        # Ns where q<mask_qual, or all Ns if mean q<min_qual, X if missing
        p7_seq = ""

        if settings["dualindex"]:
            if hasattr(res, "p7_seq") and hasattr(res, "p7_qual"):
                p7_seq = (
                    res.p7_seq
                    if res.p7_qual >= settings["mean_qual"]
                    else "N" * settings["tlen_p7"]
                )
            else:
                p7_seq = "X" * settings["tlen_p7"]

        idx_seq = p5_seq + p7_seq

        if idx_seq not in res_pfq:
            res_pfq[idx_seq] = {}

        if res.bc_seq not in res_pfq[idx_seq]:
            res_pfq[idx_seq][res.bc_seq] = 0

        res_pfq[idx_seq][res.bc_seq] += 1

    return res_pfq


def stream_fastq(file_path: str) -> Generator[Tuple[str, str, str, str], None, None]:
    """
    Iteratively stream fq or fq.gz entries into memory as 4-line string.
    Shows an indeterminate tqdm progress bar (no pre-counting).
    """
    gz = file_path.endswith(".gz")

    if gz:
        # Buffer the underlying fileobj; gzip itself then reads bigger chunks
        raw = open(file_path, "rb", buffering=1024 * 1024)
        fh = gzip.open(raw, "rt", newline="")
    else:
        fh = open(file_path, "rt", buffering=1024 * 1024, newline="")

    with (
        fh as fastq_file,
        tqdm(desc="Counting index/barcode pairs", unit="rec", ascii="-=") as pbar,
    ):
        it = iter(fastq_file)
        update_every = 1024
        n = 0

        while True:
            h = next(it, None)
            if h is None:
                break
            s = next(it)
            p = next(it)
            q = next(it)

            # cheaper than .strip(); we only want to drop newline
            yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")

            n += 1
            if n % update_every == 0:
                pbar.update(update_every)

        # flush remainder
        rem = n % update_every
        if rem:
            pbar.update(rem)


def mean_q(quals: list[int]) -> float:
    """
    Tiny hot-path helper: avoid NumPy overhead on short per-read quality slices.
    """
    return 0.0 if not quals else sum(quals) / len(quals)


class parse_fastq_rec(object):
    """
    Class to handle each fastq record as an instance while it streams.
    Return per-read idx (named if expected), bc, and assoc mean qualities.
    """

    __slots__ = (
        "seq",
        "qual",
        "p5_seq",
        "p5_qual",
        "bc_seq",
        "bc_qual",
        "p7_seq",
        "p7_qual",
    )

    # Class-level constants — built once, shared globally
    _RC_TABLE = None
    _PAT_P5 = None
    _PAT_P7 = None
    _PAT_BC = None

    @classmethod
    def init_patterns(cls, settings: dict):
        """
        Call once before streaming to precompile fuzzy patterns.
        This is much faster.
        """

        # Get allowed mismatches
        mm = settings["mismatches"]

        # compile regex patterns once; BESTMATCH yields the minimal-edit match
        cls._PAT_P5 = regex.compile(
            rf"(?:{regex.escape(settings['ref_p5'])}){{s<={mm}}}"  # FYI, e<={mm} & flags=regex.BESTMATCH takes somewhat longer but may be cleaner...
        )

        cls._PAT_BC = regex.compile(
            rf"(?:{regex.escape(settings['ref_bc'])}){{s<={mm}}}"
        )

        if settings.get("dualindex"):
            cls._PAT_P7 = regex.compile(
                rf"(?:{regex.escape(settings['ref_p7'])}){{s<={mm}}}"
            )

    def __init__(
        self,
        seq: str,
        qual: str,
        settings: dict,
    ):

        # Initialize sequence and quality
        self.seq = seq
        self.qual = self.Q33conv(qual)

        # (Opt) mask lowQ bases
        if settings["mask"]:
            self.seq = "".join(
                "N" if q < settings["mask_qual"] else b
                for b, q in zip(self.seq, self.qual)
            )

        # Process the fastq record according to settings
        self.process_sequence(settings)

    def process_sequence(self, settings):
        """
        Extract regions of interest and assoc. quality strings and update class attrs.
        Final attrs:
        - seq, qual
        - p5_seq, p5_qual
        - bc_seq, bc_qual
        - (Opt) p7_seq, p7_qual
        """
        # Extract p5 index
        p5 = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            pat=cast("regex.Pattern", parse_fastq_rec._PAT_P5),
            targ_dir=settings["tdir_p5"],
            targ_len=settings["tlen_p5"],
        )

        # If p5 found, set attrs and continue
        if p5[0]:
            self.p5_seq, self.p5_qual = p5[0], mean_q(p5[1])

            # Reverse complement if called
            if settings["rdir"] == "rev":
                self.seq = self.revcomp(self.seq)
                self.qual = self.qual[::-1]

            # Process as single or dual index
            if not settings["dualindex"]:
                self.process_single_index(settings)
            else:
                self.process_dual_index(settings)

    def process_single_index(self, settings):
        """
        Obtain barcode based on target and length
        """
        # Extract barcode for single index
        bc = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            pat=cast("regex.Pattern", parse_fastq_rec._PAT_BC),
            targ_dir=settings["tdir_bc"],
            targ_len=settings["tlen_bc"],
        )

        if bc[0]:
            self.bc_seq, self.bc_qual = bc[0], mean_q(bc[1])

    def process_dual_index(self, settings):
        """
        Obtain p7, and if present obtain barcode based on two targets to
        facilitate optional repair of short barcode, given tlen_bc and a fill sequence.
        """
        # Get di_profile
        base_profile = settings.get("base_profile")

        # Extract p7 index
        p7 = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            pat=cast("regex.Pattern", parse_fastq_rec._PAT_P7),
            targ_dir=settings["tdir_p7"],
            targ_len=settings["tlen_p7"],
        )

        # Exit loop if no p7
        if not p7[0]:
            return

        # If found, set attrs and continue
        self.p7_seq, self.p7_qual = p7[0], mean_q(p7[1])

        # Finish based on di_profile
        # Extract barcode for dual index:

        if base_profile == "X":
            # `X` has a long amplicon, just extract downstream from (pre-compiled) ref_bc. Also...
            bc = self.match_extract_single_targ(
                seq=self.seq,
                qual=self.qual,
                pat=cast("regex.Pattern", parse_fastq_rec._PAT_BC),
                targ_dir=settings["tdir_bc"],
                targ_len=settings["tlen_bc"],
            )
            # The X/INT p7 seq is revcomp, but bc is forward. I don't like this fix, but it works.
            self.p7_seq = self.revcomp(self.p7_seq)

        else:
            bc = self.match_extract_dual_targ(  # `M`s are flanked by PAT_P7 and PAT_BC, cut between
                seq=self.seq,
                qual=self.qual,
                pat1=cast("regex.Pattern", parse_fastq_rec._PAT_P7),
                pat2=cast("regex.Pattern", parse_fastq_rec._PAT_BC),
                targ_len=settings["tlen_bc"],
                fill_seq=settings["fill_seq"],
            )

        if bc[0]:
            self.bc_seq, self.bc_qual = bc[0], mean_q(bc[1])

    @staticmethod
    def Q33conv(qual: str):
        """
        Convert a string of Q33 ascii characters into a list of PHRED scores.
        """
        return [ord(q) - 33 for q in qual]

    @staticmethod
    def revcomp(seq: str):
        """
        Returns the reverse complement of a DNA sequence.
        Lazily initialize the RC table if needed.
        """
        if parse_fastq_rec._RC_TABLE is None:
            parse_fastq_rec._RC_TABLE = str.maketrans(
                {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "X": "X"}
            )
        return seq.translate(parse_fastq_rec._RC_TABLE)[::-1]

    @staticmethod
    def match_extract_single_targ(
        seq: str,
        qual: list[int],
        pat: "regex.Pattern",
        targ_dir: str,
        targ_len: int,
    ):
        """
        Extract sequence of given length and direction by single target.
        """
        # locate match
        m = pat.search(seq)

        if m:
            # get positions
            start, end = m.start(), m.end()

            # return appropriate slice
            if targ_dir == "upstream":
                return [seq[start - targ_len : start], qual[start - targ_len : start]]
            elif targ_dir == "downstream":
                return [seq[end : end + targ_len], qual[end : end + targ_len]]

            # or tell user we need a direction
            else:
                raise ValueError("targ_dir param requires 'upstream' or 'downstream'")

        # if no match, return empty
        return ["", []]

    @staticmethod
    def match_extract_dual_targ(
        seq: str,
        qual: list[int],
        pat1: "regex.Pattern",
        pat2: "regex.Pattern",
        targ_len: int,
        fill_seq: str | None,
    ):
        """
        Extract sequence between ref1 and ref2.
        Pad front of xseq to targ_len with 3' of fill_seq, if fill_seq provided
        """
        # locate matches
        m1 = pat1.search(seq)
        m2 = pat2.search(seq)

        if m1 and m2:
            start, end = m1.end(), m2.start()

            # if start is before or at end (this catches empty short bcs!)
            if start <= end:
                # get slices
                xseq, xqual = seq[start:end], qual[start:end]

                # pad 5' bc_seq to desired length with 3' of fill_seq
                # (if fill_seq is provided in stock settings, and bc is short)
                if fill_seq and len(xseq) < targ_len:
                    diff = targ_len - len(xseq)
                    xseq = fill_seq[-diff:] + xseq
                    xqual = [30] * diff + xqual

                return [xseq, xqual]

        return ["", []]


# %% versions
"""
v1.0.0 20260625
    - Created as the linear parse-mode pipeline module.
    - Combines former ingest/countdict.py reference-loading, count-dict,
      dataframe, runinfo-merge, QC, and SeqRun construction machinery.
    - Combines former ingest/fastq.py FASTQ streaming and per-record parsing
      machinery.
    - Keeps workbook/csv emission outside ParseMode in __main__/emit.
"""
