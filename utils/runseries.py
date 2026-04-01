# runseries.py

"""
Name:       runseries.py
Author:     CAG
Version:    1.0
Date:       2026/03/26

`RunSeries` is the top-level model used by compile mode:
- n `SeqRun` objects + a cached raw compiled long dataframe
- helpers for ingest, compile, QC, and emission
"""

#%% Imports

from __future__ import annotations

from dataclasses import dataclass, field
import os
from typing import TYPE_CHECKING, Any, Iterable, Mapping

import pandas as pd

import utils.qc_ops as qc
from utils.runinfo import RunInfo
from utils.seqrun import SeqRun

if TYPE_CHECKING:
    from utils.parse_xlsx import ParsedAnalysis

#%% RunSeries model

@dataclass
class RunSeries:
    """
    Container for a series of sequencing runs.
    Handles multi-run compilation, cross-run QC, and output emission for compile mode.

    Identity contract for compile mode:
    - within a run, samples are identified by `idx_name` when available
    - across runs, callers should distinguish samples by `run_number + idx_name`
    - legacy/base-csv rows without idx_name fall back to the composite
      group/name/date/run_number identity
    """
    runs: dict[str, SeqRun] = field(default_factory=dict)      # Primary structured representation: one SeqRun per run identity.
    series_meta: dict[str, Any] = field(default_factory=dict)  # Top-level build metadata for the series.
    data_df: pd.DataFrame | None = None                        # Cached flattened raw long series df.
    runinfo_by_run_id: dict[str, RunInfo] = field(default_factory=dict)
    source: str = "unknown"

    def __post_init__(self) -> None:
        # Make private copies of mutable inputs to protect source
        self.runs = self._coerce_runs(self.runs)
        self.series_meta = dict(self.series_meta)
        self.data_df = None if self.data_df is None else self.data_df.copy()
        self.runinfo_by_run_id = {
            str(run_id): runinfo
            for run_id, runinfo in self.runinfo_by_run_id.items()
            if runinfo is not None
        }

    # ====================
    # INGEST
    # ====================

    @classmethod
    def from_long_df(
        cls,
        df: pd.DataFrame,
        *,
        source: str = "unknown",
        runinfo_by_run_id: Mapping[str, RunInfo] | None = None,
    ) -> "RunSeries":
        """
        Materialize a RunSeries from a long-format dataframe spanning one or more runs.

        This is the main dataframe-to-container ingest path for compile mode.
        The input df is assumed to already be normalized and to represent raw
        compiled series data. This method's job is to split that long df into
        per-run slices and rebuild one `SeqRun` for each run identity.
        """
        df = qc.normalize_long_df(df)

        runinfo_by_run_id = {} if runinfo_by_run_id is None else dict(runinfo_by_run_id)

        if df.empty:
            return cls(
                runinfo_by_run_id=runinfo_by_run_id,
                source=source,
            )

        # Choose grouping column for run-level slicing.
        group_col = SeqRun.long_df_run_group_col(df)

        if group_col is None:
            # No usable run-level grouping column was found, so treat the whole
            # dataframe as one run and key it with a fallback id.
            return cls(
                runs={
                    "0": SeqRun.from_long_df(
                        df,
                        source=source,
                        runinfo=runinfo_by_run_id.get("0"),
                    )
                },
                series_meta={"groupby_col": None},
                runinfo_by_run_id=runinfo_by_run_id,
                source=source,
            )

        runs = {
            ("0" if pd.isna(run_id) else str(run_id)): SeqRun.from_long_df(
                run_df,
                source=source,
                runinfo=runinfo_by_run_id.get("0" if pd.isna(run_id) else str(run_id)),
            )
            for run_id, run_df in df.groupby(group_col, dropna=False, sort=False)
        }

        return cls(
            runs=runs,
            series_meta={"groupby_col": group_col},
            runinfo_by_run_id=runinfo_by_run_id,
            source=source,
        )


    @staticmethod
    def _coerce_runs(runs: Mapping[str, SeqRun] | Iterable[SeqRun]) -> dict[str, SeqRun]:
        """
        Normalize incoming run containers to the internal `{run_id: SeqRun}` shape.

        This is used during `RunSeries` initialization so callers can pass either:
        - a dict that is already keyed by run id, or
        - any iterable of `SeqRun` objects

        If no explicit key is provided, a fallback run id is derived from the run.
        """
        if isinstance(runs, Mapping):
            items = runs.items()
        else:
            items = ((None, run) for run in runs)

        run_map: dict[str, SeqRun] = {}
        for key, run in items:
            # Prefer an explicit mapping key. Otherwise derive a stable-ish id
            # from run metadata so the series still has addressable runs.
            run_id = key or str(run.run_number or run.run_name or run.filename or len(run_map))
            if run_id in run_map:
                raise ValueError(f"Duplicate run_id in RunSeries: {run_id}")
            run_map[str(run_id)] = run
        return run_map


    # ====================
    # COMPILE
    # ====================

    @staticmethod
    def read_base_csv(base_csv_path: str) -> pd.DataFrame:
        """
        Read and lightly normalize a previously compiled all.csv file.
        """
        return (
            pd.read_csv(base_csv_path)
            .drop(columns=["ldist_group_lvl", "ldist_all"], errors="ignore")
        )

    @classmethod
    def build_from_parsed_files(
        cls,
        *,
        parsed_files: list[pd.DataFrame | SeqRun | ParsedAnalysis],
        settings: dict,
        base_df: pd.DataFrame | None = None,
    ) -> "RunSeries":
        """
        Build one raw deduplicated RunSeries from parsed analyses and optional existing compiled data.
        """
        # Normalize and deduplicate compile inputs.
        df, runinfo_by_run_id, series_meta = cls._compile_inputs(
            parsed_files=parsed_files,
            existing_df=base_df,
            out_prefix=settings.get("out_prefix"),
        )

        # Materialize the raw compiled series.
        run_series = cls.from_long_df(
            df,
            source="compiled",
            runinfo_by_run_id=runinfo_by_run_id,
        )

        run_series.series_meta.update(series_meta)
        run_series.data_df = df.copy()

        return run_series


    @classmethod
    def _compile_inputs(
        cls,
        *,
        parsed_files: list[pd.DataFrame | SeqRun | ParsedAnalysis],
        existing_df: pd.DataFrame | None = None,
        out_prefix: str | None = None,
    ) -> tuple[pd.DataFrame, dict[str, RunInfo], dict[str, Any]]:
        """
        Normalize compile inputs and resolve them to one deduplicated long dataframe.
        """
        parsed_df, runinfo_by_run_id = cls._coerce_parsed_inputs(parsed_files)
        existing_df = cls._normalize_existing_df(existing_df)

        print("Compile: deduplicating parsed analyses and existing data...")
        combined_df = cls._deduplicate_compile_sources(
            parsed_df=parsed_df,
            existing_df=existing_df,
        )

        print("Compile: assembling final combined dataframe...")
        combined_df = cls._sort_df(combined_df)
        combined_df = cls._disambiguate_sample_ids(combined_df)

        print("Compile: running core QC...")
        combined_df = cls._apply_core_qc(combined_df)

        n_source_runs = 0
        if not combined_df.empty:
            run_col = "run_id" if "run_id" in combined_df.columns else "run_number"
            n_source_runs = combined_df[run_col].dropna().astype(str).nunique()

        return (
            combined_df,
            runinfo_by_run_id,
            {
                "out_prefix": out_prefix,
                "n_source_runs": n_source_runs,
            },
        )


    @staticmethod
    def _coerce_parsed_inputs(
        parsed_files: list[pd.DataFrame | SeqRun | ParsedAnalysis],
    ) -> tuple[pd.DataFrame, dict[str, RunInfo]]:
        """
        Normalize compile-mode inputs to one long dataframe plus a runinfo lookup.

        This is the first compile step used by `_compile_inputs()`.
        It accepts mixed input types from compile mode, converts each to the
        normalized long-df shape, records runinfo by run_id when available, and
        stamps parsed analysis rows with source metadata for later deduplication.
        """
        parsed_dfs: list[pd.DataFrame] = []
        runinfo_by_run_id: dict[str, RunInfo] = {}

        for run_order, item in enumerate(parsed_files):
            runinfo = None
            filepath = None

            # Accept parsed workbook objects, SeqRun objects, or bare dataframes.
            if isinstance(item, SeqRun):
                df = item.df.copy()
                runinfo = item.runinfo
                filepath = item.filepath
            elif hasattr(item, "data_df") and hasattr(item, "runinfo") and hasattr(item, "filepath"):
                if item.data_df is None:
                    continue
                df = item.data_df.copy()
                runinfo = item.runinfo
                filepath = item.filepath
            else:
                df = item.copy()

            df = qc.normalize_long_df(df)
            if df.empty:
                continue

            run_ids = df["run_id"].dropna().astype(str).unique().tolist()
            if not run_ids:
                raise ValueError("Parsed compile input is missing run_id.")
            if len(run_ids) != 1:
                raise ValueError(f"Expected one run_id per parsed analysis, found {run_ids}")

            run_id = run_ids[0]

            # Stamp analysis rows so dedup can prefer newer parsed inputs over csv rows.
            parsed_dfs.append(
                RunSeries._stamp_compile_source(
                    df=df,
                    source_label="analysis",
                    source_rank=1,
                    source_order=run_order,
                    source_mtime=RunSeries._filepath_mtime(filepath),
                )
            )

            if runinfo is not None:
                runinfo_by_run_id[run_id] = runinfo

        if not parsed_dfs:
            return pd.DataFrame(), runinfo_by_run_id

        return pd.concat(parsed_dfs, ignore_index=True), runinfo_by_run_id


    @staticmethod
    def _normalize_existing_df(existing_df: pd.DataFrame | None) -> pd.DataFrame:
        """
        Normalize an optional existing compiled long dataframe.
        """
        if existing_df is None:
            return pd.DataFrame()
        return qc.normalize_long_df(existing_df)


    @staticmethod
    def _filepath_mtime(filepath: Any) -> float:
        """
        Best-effort mtime lookup used for parsed-analysis precedence.
        """
        if filepath is None:
            return float("-inf")

        try:
            return filepath.stat().st_mtime
        except (AttributeError, OSError):
            return float("-inf")


    @staticmethod
    def _stamp_compile_source(
        *,
        df: pd.DataFrame,
        source_label: str,
        source_rank: int,
        source_order: int,
        source_mtime: float,
    ) -> pd.DataFrame:
        """
        Attach temporary source-priority columns used during compile deduplication.
        """
        stamped = qc.normalize_long_df(df).copy()
        stamped["__compile_source"] = source_label
        stamped["__compile_source_rank"] = source_rank
        stamped["__compile_source_order"] = source_order
        stamped["__compile_source_mtime"] = source_mtime
        return stamped

    @staticmethod
    def _describe_compile_sample(row: pd.Series) -> str:
        """
        Build a readable sample label for compile dedup/drop messages.
        """
        parts: list[str] = []

        run_number = row.get("run_number")
        if pd.notna(run_number) and str(run_number).strip():
            parts.append(f"run {run_number}")

        idx_name = row.get("idx_name")
        if pd.notna(idx_name) and str(idx_name).strip():
            parts.append(f"idx {idx_name}")

        sample_id = row.get("sample_id")
        if pd.notna(sample_id) and str(sample_id).strip():
            parts.append(f"sample_id {sample_id}")

        samp_group = row.get("samp_group")
        samp_name = row.get("samp_name")
        samp_date = row.get("samp_date")
        sample_bits = [
            str(value).strip()
            for value in [samp_group, samp_name]
            if pd.notna(value) and str(value).strip()
        ]
        if pd.notna(samp_date):
            norm_date = qc.normalize_date(samp_date)
            if pd.notna(norm_date):
                sample_bits.append(str(norm_date.date()))
            elif str(samp_date).strip():
                sample_bits.append(str(samp_date).strip())
        if sample_bits:
            parts.append(" / ".join(sample_bits))

        source_label = row.get("__compile_source")
        if pd.notna(source_label) and str(source_label).strip():
            parts.append(f"source={source_label}")

        filename = row.get("filename")
        if pd.notna(filename) and str(filename).strip():
            parts.append(f"file={filename}")

        return " | ".join(parts) if parts else "<unlabeled sample>"

    @classmethod
    def _report_dropped_compile_replicates(
        cls,
        *,
        dedup_df: pd.DataFrame,
        winners: pd.DataFrame,
    ) -> None:
        """
        Print explicit compile dedup decisions so dropped replicates are visible.
        """
        if dedup_df.empty:
            return

        winner_rows = winners.rename(
            columns={
                "__compile_source_rank": "__winner_source_rank",
                "__compile_source_mtime": "__winner_source_mtime",
                "__compile_source_order": "__winner_source_order",
            }
        )

        decisions = dedup_df.merge(
            winner_rows,
            how="left",
            on=["__sample_identity"],
        )

        kept_mask = (
            decisions["__compile_source_rank"].eq(decisions["__winner_source_rank"])
            & decisions["__compile_source_mtime"].eq(decisions["__winner_source_mtime"])
            & decisions["__compile_source_order"].eq(decisions["__winner_source_order"])
        )
        dropped = decisions.loc[~kept_mask].copy()
        if dropped.empty:
            return

        print("Compile: dropping duplicate replicate inputs...")
        for sample_identity, group in dropped.groupby("__sample_identity", sort=False):
            kept_row = decisions.loc[decisions["__sample_identity"].eq(sample_identity) & kept_mask]
            kept_desc = cls._describe_compile_sample(kept_row.iloc[0]) if not kept_row.empty else "<unknown kept source>"
            print(f"  keeping: {kept_desc}")
            for _, row in group.iterrows():
                print(f"  dropped: {cls._describe_compile_sample(row)}")


    @classmethod
    def _deduplicate_compile_sources(
        cls,
        *,
        parsed_df: pd.DataFrame,
        existing_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Deduplicate compile inputs before QC using sample identity.

        Prefer `run_number + idx_name` when idx_name exists. Fall back to the
        older metadata composite only when idx_name is unavailable, which is
        expected mainly for existing/base csv rows.
        """
        stamped_frames: list[pd.DataFrame] = []

        # Existing compiled csv participates in dedup with lower precedence.
        if not existing_df.empty:
            stamped_frames.append(
                cls._stamp_compile_source(
                    df=existing_df,
                    source_label="csv",
                    source_rank=0,
                    source_order=-1,
                    source_mtime=float("-inf"),
                )
            )

        if not parsed_df.empty:
            stamped_frames.append(parsed_df.copy())

        if not stamped_frames:
            return pd.DataFrame()

        combined = pd.concat(stamped_frames, ignore_index=True)
        combined = qc.normalize_long_df(combined)

        required_cols = ["run_number", "idx_name", "samp_group", "samp_name", "samp_date"]
        missing_cols = [col for col in required_cols if col not in combined.columns]
        if missing_cols:
            raise ValueError(
                f"Compile dedup requires columns {required_cols}; missing {missing_cols}"
            )

        dedup_df = combined.copy()
        dedup_df["run_number"] = dedup_df["run_number"].astype(str)
        dedup_df["idx_name"] = dedup_df["idx_name"].astype("string")
        dedup_df["samp_group"] = dedup_df["samp_group"].astype(str)
        dedup_df["samp_name"] = dedup_df["samp_name"].astype(str)
        dedup_df["samp_date"] = dedup_df["samp_date"].apply(qc.normalize_date)
        dedup_df["__sample_identity"] = dedup_df.apply(
            lambda row: (
                f"{row.get('run_number')}::{row.get('idx_name')}"
                if pd.notna(row.get("idx_name")) and str(row.get("idx_name")).strip()
                else (
                    f"{row.get('run_number')}::{row.get('sample_id')}"
                    if pd.notna(row.get("sample_id")) and str(row.get("sample_id")).strip()
                    else f"{row.get('run_number')}::{row.get('samp_group')}::{row.get('samp_name')}::{row.get('samp_date')}"
                )
            ),
            axis=1,
        )

        # One winning source per logical sample identity. Parsed analyses beat
        # existing/base csv rows when both describe the same sample.
        winners = (
            dedup_df[
                ["__sample_identity"]
                + [
                    "__compile_source_rank",
                    "__compile_source_mtime",
                    "__compile_source_order",
                ]
            ]
            .drop_duplicates()
            .sort_values(
                by=["__sample_identity"]
                + [
                    "__compile_source_rank",
                    "__compile_source_mtime",
                    "__compile_source_order",
                ],
                ascending=[True, True, True, True],
                kind="stable",
            )
            .drop_duplicates(subset=["__sample_identity"], keep="last")
        )

        cls._report_dropped_compile_replicates(
            dedup_df=dedup_df,
            winners=winners,
        )

        merged = dedup_df.merge(
            winners,
            how="inner",
            on=["__sample_identity"]
            + [
                "__compile_source_rank",
                "__compile_source_mtime",
                "__compile_source_order",
            ],
        )

        return merged.drop(
            columns=[
                "__compile_source",
                "__compile_source_rank",
                "__compile_source_mtime",
                "__compile_source_order",
                "__sample_identity",
            ],
            errors="ignore",
        )

    @staticmethod
    def _disambiguate_sample_ids(df: pd.DataFrame) -> pd.DataFrame:
        """
        Disambiguate repeated compile-mode sample_id values as ::dupN by first encounter.

        Distinct sample blocks are identified using sample-level metadata, including
        idx_name when present, but the emitted suffix is always dupN.
        """
        if df.empty or "sample_id" not in df.columns:
            return df

        df = df.copy()
        df["sample_id"] = df["sample_id"].fillna("").astype(str)

        sample_level_cols = [
            "sample_id",
            "idx_name",
            "run_id",
            "run_number",
            "run_name",
            "run_date",
            "samp_group",
            "samp_name",
            "samp_date",
            "input",
            "filename",
        ]
        present_cols = [col for col in sample_level_cols if col in df.columns]
        if not present_cols:
            return df

        sample_blocks = df.loc[:, present_cols].drop_duplicates().reset_index(drop=True)
        if sample_blocks.empty:
            return df

        sample_blocks["__dup_rank"] = sample_blocks.groupby("sample_id", sort=False).cumcount() + 1
        has_dups = sample_blocks.groupby("sample_id", sort=False)["sample_id"].transform("size") > 1
        sample_blocks["__resolved_sample_id"] = sample_blocks["sample_id"]
        sample_blocks.loc[has_dups, "__resolved_sample_id"] = (
            sample_blocks.loc[has_dups, "sample_id"]
            + "::dup"
            + sample_blocks.loc[has_dups, "__dup_rank"].astype(str)
        )

        dup_blocks = sample_blocks.loc[has_dups, ["sample_id", "__resolved_sample_id"]].copy()
        if not dup_blocks.empty:
            print("Compile: disambiguating duplicated sample_id values...")
            for sample_id, rows in dup_blocks.groupby("sample_id", sort=False):
                resolved_ids = rows["__resolved_sample_id"].tolist()
                print(f"  {sample_id} -> {resolved_ids}")

        merged = df.merge(
            sample_blocks[present_cols + ["__resolved_sample_id"]],
            how="left",
            on=present_cols,
        )
        resolved_mask = merged["__resolved_sample_id"].notna()
        merged.loc[resolved_mask, "sample_id"] = merged.loc[resolved_mask, "__resolved_sample_id"]
        return merged.drop(columns=["__resolved_sample_id"])


    @staticmethod
    def _sort_df(df: pd.DataFrame) -> pd.DataFrame:
        """
        Prioritize key metadata columns and sort rows by sample identity.
        """
        df = df.copy()
        cols = ["samp_group", "samp_name", "samp_date"]

        df = df[[c for c in cols if c in df.columns] + [c for c in df.columns if c not in cols]]
        df = df.sort_values(
            by=cols,
            ascending=[True] * len(cols),
            kind="stable",
        ).reset_index(drop=True)

        return df


    @staticmethod
    def _apply_core_qc(df: pd.DataFrame) -> pd.DataFrame:
        """
        Run core QC operations on the combined compile dataframe.
        """
        df = df.copy()

        cols = ["bc_count", "proportion", "input"]
        df[cols] = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0)

        df = qc.rank_uniques(df=df)
        df["bc_name"] = df["bc_name"].str.replace("SIVmac293M2", "SIVmac239M2", regex=False)
        df = qc.flag_shorts(df=df)

        df["ldist_samp_lvl"] = df["ldist_samp_lvl"].str.replace(r"^(1-)|(-1)$", "", regex=True)
        df["ldist_samp_lvl"] = df["ldist_samp_lvl"].str.replace(r"^(Unique)", "", regex=True)
        df["ldist_samp_lvl"] = df["ldist_samp_lvl"].replace("NA-NA", None)

        df["multi_idx"] = (
            df.groupby(["run_number", "samp_group", "bc_name"])["samp_name"]
            .transform(lambda s: s.nunique() > 1)
        )

        df = SeqRun.annotate_x_idx(df)

        return df

    # ====================
    # EMIT
    # ====================

    @property
    def run_list(self) -> list[SeqRun]:
        return list(self.runs.values())

    @property
    def run_ids(self) -> list[str]:
        return list(self.runs.keys())

    @property
    def n_runs(self) -> int:
        return len(self.runs)

    @property
    def df(self) -> pd.DataFrame:
        """
        Flatten all runs into the raw compiled long-format dataframe.
        """
        if self.data_df is not None:
            return self.data_df.copy()
        if not self.runs:
            return pd.DataFrame()
        return pd.concat([run.df for run in self.run_list], ignore_index=True)


    def copy(self, **updates: Any) -> "RunSeries":
        """
        Deep-ish copy used when callers want another mutable series object.
        """
        data = {
            "runs": {run_id: run.copy() for run_id, run in self.runs.items()},
            "series_meta": self.series_meta.copy(),
            "data_df": None if self.data_df is None else self.data_df.copy(),
            "runinfo_by_run_id": {
                run_id: RunInfo(
                    run_meta=runinfo.run_meta.copy(),
                    sample_meta={k: v.copy() for k, v in runinfo.sample_meta.items()},
                    raw_df=None if runinfo.raw_df is None else runinfo.raw_df.copy(),
                    source=runinfo.source,
                    filepath=runinfo.filepath,
                )
                for run_id, runinfo in self.runinfo_by_run_id.items()
            },
            "source": self.source,
        }
        data.update(updates)
        return RunSeries(**data)


    def write_output(self, settings: dict) -> None:
        """
        Format and save compiled series output files.
        """
        from utils.data_views import (
            concat_series_runinfo,
            contam_report,
            group_mat_dict,
            sidelong_tables,
        )

        # Set up paths and filenames
        out_path   = settings["out_path"]
        out_prefix = settings["out_prefix"]

        os.makedirs(out_path, exist_ok=True)
        csv_path = os.path.join(out_path, f"{out_prefix}_all.csv")
        xlsx_path = os.path.join(out_path, f"{out_prefix}_above_cutoff.xlsx")

        # Build QC-aware compile outputs before writing files so all.csv
        # includes parent-check annotations produced during emit_qc().
        compile_df, df_ac, contam_check = self.emit_qc(settings=settings)
        compile_df_out = qc.format_output_dates(compile_df, columns=["run_date", "samp_date"])
        df_ac_out = qc.format_output_dates(df_ac, columns=["run_date", "samp_date"])

        # Write full compiled df to csv
        compile_df_out.to_csv(csv_path, index=False)

        # Format workbook components
        runinfo_df = qc.format_output_dates(
            concat_series_runinfo(self),
            columns=["Date", "run_date", "samp_date"],
        )
        sidelong_dict, _ = sidelong_tables(df_ac_out)
        _, mat_dict = group_mat_dict(df_ac_out)

        # Write the above-cutoff workbook
        with pd.ExcelWriter(xlsx_path) as writer:
            if not runinfo_df.empty:
                runinfo_df.to_excel(writer, sheet_name="runinfo", index=False)

            # above_cutoff 
            df_ac_out.to_excel(writer, sheet_name="all_samples", index=False)

            # full matrix and contamination report
            if contam_check:
                full_mat, contam_table = contam_report(df_ac_out)
                full_mat.to_excel(writer, sheet_name="full_matrix", index=True, header=True)
                contam_table.to_excel(writer, sheet_name="contam_report", index=True, header=True)

            # matrix and sidelongs
            for samp_group in mat_dict:
                print(f"Writing {samp_group} ...")
                sidelong_dict[samp_group].to_excel(
                    writer, sheet_name=f"{samp_group} samples", index=False, header=False
                )
                mat_dict[samp_group].to_excel(
                    writer, sheet_name=f"{samp_group} matrix", index=True, header=True
                )


    def emit_qc(
        self,
        *,
        settings: dict,
    ) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
        """
        Emit QC dataframe outputs from the raw compiled series.
        """
        # Work from a copy so raw compiled state stays unchanged.
        working_series = self.copy()

        # Collapse to parent at the contained SeqSamp boundary when requested.
        if settings.get("collapse_to_parent", True):
            print("Compile: collapsing children to parent within each sample...")
            working_series.collapse_samples_in_place()

        compile_df = working_series.df.copy()
        if compile_df.empty:
            return pd.DataFrame(), pd.DataFrame(), False

        if any(settings.get(key, False) for key in ["ldist_samp_lvl", "ldist_group_lvl", "contam_check"]):
            print("Compile: running optional QC...")
        compile_df, contam_check = self._apply_optional_qc(
            df=compile_df,
            settings=settings,
        )

        print("Compile: building above-cutoff dataframe...")
        df_ac = self._build_above_cutoff_df(compile_df)
        return compile_df, df_ac, contam_check


    def collapse_samples_in_place(self) -> None:
        """
        Collapse each SeqSamp in the series to parent in place.
        """
        # Replace each sample with its collapsed version and clear cached dfs.
        for run in self.run_list:
            run.samples = {
                sample_id: sample.collapse_to_parent()
                for sample_id, sample in run.samples.items()
            }
            run.data_df = None

        self.data_df = None

    # ====================
    # QC HELPERS
    # ====================

    @classmethod
    def _apply_optional_qc(cls, df: pd.DataFrame, settings: dict) -> tuple[pd.DataFrame, bool]:
        """
        Run optional compile QC checks controlled by settings.
        """
        if settings.get("ldist_samp_lvl", False):
            df = cls._apply_ldist_qc(
                df=df,
                groupby_cols=["samp_group", "samp_name", "samp_date"],
                output_col="ldist_samp_lvl",
                drop_existing=True,
                settings=settings,
            )

        if settings.get("ldist_group_lvl", False):
            df = cls._apply_ldist_qc(
                df=df,
                groupby_cols=["samp_group"],
                aggr=True,
                output_col="ldist_group_lvl",
                settings=settings,
            )

        contam_check = bool(settings.get("contam_check", False))
        if contam_check:
            df = cls._apply_ldist_qc(
                df=df,
                groupby_cols=None,
                aggr=True,
                output_col="ldist_all",
                settings=settings,
            )

        return df, contam_check


    @staticmethod
    def _apply_ldist_qc(
        df: pd.DataFrame,
        settings: dict,
        groupby_cols: list[str] | None,
        aggr: bool = False,
        output_col: str = "putative_parent",
        drop_existing: bool = False,
    ) -> pd.DataFrame:
        """
        Run ldist QC checks with flexible grouping and output columns.
        """
        if drop_existing and output_col in df.columns:
            df = df.drop(columns=[output_col])

        if groupby_cols:
            return (
                df.groupby(groupby_cols, as_index=False)
                .apply(
                    lambda group: qc.flag_putative_parents(
                        set_name=group.name,
                        df=group,
                        settings=settings,
                        aggr=aggr,
                    )
                )
                .rename(columns={"putative_parent": output_col})
            )

        return qc.flag_putative_parents(
            set_name="Full Series",
            df=df,
            settings=settings,
            aggr=aggr,
        ).rename(columns={"putative_parent": output_col})


    @staticmethod
    def _build_above_cutoff_df(df: pd.DataFrame) -> pd.DataFrame:
        """
        Build an above-cutoff view and recalculate proportions within that scope.
        """
        if df.empty or "above_cutoff" not in df.columns:
            return df.iloc[0:0].copy()

        df_ac = df.loc[df["above_cutoff"]].copy()

        if df_ac.empty:
            return df_ac

        df_ac["bc_count"] = pd.to_numeric(df_ac["bc_count"], errors="coerce").fillna(0)
        sample_keys = ["samp_group", "samp_name", "samp_date", "run_number", "filename"]

        df_ac["__denom"] = (
            df_ac.groupby(sample_keys, dropna=False)["bc_count"].transform("sum")
        )
        df_ac["proportion"] = (df_ac["bc_count"] / df_ac["__denom"]).fillna(0)

        return df_ac.drop(columns=["__denom"])

#%% Versions
"""
v1.0 
- Initial version, accomodating new model/container structure
- Supercedes deprecated xlsx_compiler module and associated code in bcParse.py
"""
