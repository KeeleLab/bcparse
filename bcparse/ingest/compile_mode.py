# compile_mode.py
"""
Name:       compile_mode.py
Author:     CAG
Version:    1.0.0
Date:       2026/06/25

Compile-mode pipeline for a series of Analysis/Discovery workbooks with optional
prior compiled base CSV.

This module is intentionally linear, CompileMode.__init__ describes action.

The goal is to keep compile-mode behavior discoverable in one place without
turning the light container classes into alternate pipeline entry points.

Major classes and helpers:

    CompileMode
        Owns settings validation, workbook/base-CSV ingestion, deduplication,
        compile QC, RunSeries construction, and missing RunInfo synthesis.

    XlsxPathManager / AnalysisWorkbookParser
        Discover parsed workbooks and normalize them into long-format
        compiled data plus RunInfo metadata.

    ParsedAnalysis / WorkbookParseError
        Small support objects for workbook parser output and readable failures.
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

import bcparse.qc as qc
from bcparse.containers.runinfo import RunInfo
from bcparse.containers.runseries import RunSeries
from bcparse.containers.seqrun import SeqRun
from bcparse.containers.seqsamp import SeqSamp
from bcparse.qc.identity import make_run_id, normalize_date, normalize_long_df
from bcparse.settings import CompileSettings

logger = logging.getLogger(__name__)

# %% CompileMode

class CompileMode:
    """
    Top-level compile-mode orchestration for one compiled output series.

    Call site:
        bcparse.__main__._run_compile() builds this class as:
            compile_run = CompileMode(settings=compile_settings)

    Inputs:
        CompileSettings object from CLI/user configuration:
            - parse workbook directory/path settings
            - optional base compiled CSV path
            - output prefix/path settings
            - compile/QC thresholds and behavior flags

    Outputs:
        run_series:
            Final RunSeries container with SeqRun children, runinfo_by_run_id,
            compile-level series_meta, and the final compiled dataframe.
        combined_df:
            Final deduplicated and QC-normalized long dataframe used by
            RunSeries and workbook emission.
    """

    settings: CompileSettings  # Original user/CLI compile config.
    runtime_settings: dict[str, Any]  # Runtime compile settings dict.

    xlsx_path_manager: XlsxPathManager  # Parse workbook discovery helper.
    xlsx_files: list[Path]  # Parse workbook paths selected for parsing.
    parsed_files: list[ParsedAnalysis]  # Parsed non-empty parse workbooks.
    base_df: pd.DataFrame | None  # Optional prior compiled CSV dataframe.

    compile_frames: list[pd.DataFrame]  # Source-ranked long dataframe inputs.
    runinfo_by_run_id: dict[str, RunInfo]  # RunInfo objects keyed by run_id.
    combined_df: pd.DataFrame  # Final deduplicated/QC-normalized dataframe.

    n_source_runs: int  # Number of run_id values in combined_df.
    run_series: RunSeries  # Final compile-mode RunSeries container.

    def __init__(self, *, settings: CompileSettings) -> None:

        self.settings = settings
        self.settings.validate()
        self.runtime_settings = settings.to_dict()
        self.runtime_settings["core_bc"] = (
            None  # For now, until I figure out how to detect stocks and do core search...
        )

        # 1. Discover parse workbooks.
        self.xlsx_path_manager = XlsxPathManager(
            wdir=self.runtime_settings["xlsx_path"]
        )
        self.xlsx_files = self.xlsx_path_manager.files
        print(f"Discovered {len(self.xlsx_files)} parse workbook(s).")

        # 2. Parse each workbook.
        self.parsed_files = []
        self._parse_analysis_workbooks()

        # 3. Read optional base CSV.
        self.base_df = self._read_base_csv()

        # 4. Normalize all compile dataframes and attach source provenance.
        print("Compiling series...")
        self.compile_frames, self.runinfo_by_run_id = self._build_frames()

        # 5. Deduplicate by sample_id.
        #    Current policy:
        #    - analysis beats csv
        #    - newer mtime beats older input
        #    - later run order breaks remaining ties
        print("- deduplicating parsed analyses and existing data...")
        self.combined_df = self._dedup_sample(
            pd.concat(self.compile_frames, ignore_index=True)
        )
        self.combined_df = self._assemble_final_combined_dataframe()

        # 6. Run core QC.
        print("- running core QC...")
        self.combined_df = self._core_qc()

        # 7. Build RunSeries.
        self.n_source_runs = self.combined_df["run_id"].dropna().astype(str).nunique()
        self.run_series = self._build_runseries(source="compiled")

        # 8. Synthesize RunInfo for base-CSV-derived runs.
        self._ensure_missing_runinfo()

        # 9. Store compile-level metadata.
        self._store_series_meta()


    def _parse_analysis_workbooks(self) -> None:
        """
        Parse discovered workbooks and retain non-empty parsed results.
        """
        print()
        for i, filepath in enumerate(self.xlsx_files, start=1):
            print(f"Parsing workbook {i}/{len(self.xlsx_files)}: {filepath.name}")
            try:
                parsed = AnalysisWorkbookParser(filepath=filepath)
                if not parsed.data_df.empty:
                    self.parsed_files.append(parsed.analysis)
            except WorkbookParseError as e:
                raise SystemExit(f"Compile aborted while parsing workbook:\n{e}") from e
        print()


    def _read_base_csv(self) -> pd.DataFrame | None:
        """
        Load the optional prior compiled CSV declared in runtime settings.
        """
        base_csv_path = self.runtime_settings.get("base_csv_path")
        if not base_csv_path:
            return None

        base_df = self._read_base_csv_file(str(base_csv_path))
        print(
            f"Loaded base csv: {Path(base_csv_path).name} with {len(base_df)} row(s)."
        )
        return base_df


    def _read_base_csv_file(self, base_csv_path: str) -> pd.DataFrame:
        """
        Read prior compiled CSV rows, dropping columns recomputed by this run.
        """
        return pd.read_csv(base_csv_path).drop(
            columns=["ldist_group_lvl", "ldist_all"], errors="ignore"
        )


    def _build_frames(
        self,
    ) -> tuple[list[pd.DataFrame], dict[str, RunInfo]]:
        """
        Build source-ranked long-dataframe inputs from base CSV and analyses.
        """
        compile_frames: list[pd.DataFrame] = []
        runinfo_by_run_id: dict[str, RunInfo] = {}

        if self.base_df is not None:
            csv_df = qc.normalize_long_df(
                self.base_df.drop(columns=["sample_id"], errors="ignore"),
                source_path=self.runtime_settings.get("base_csv_path"),
            )
            if not csv_df.empty:
                csv_df["__compile_source"] = "csv"
                csv_df["__compile_source_note"] = self._source_note(
                    "base_csv",
                    self.runtime_settings.get("base_csv_path"),
                )
                csv_df["__compile_source_rank"] = 0
                csv_df["__compile_source_order"] = -1
                csv_df["__compile_source_mtime"] = float("-inf")
                compile_frames.append(csv_df)

        for run_order, parsed in enumerate(self.parsed_files):
            df = parsed.data_df.copy()
            if df.empty:
                continue

            run_ids = df["run_id"].dropna().astype(str).unique().tolist()
            if len(run_ids) != 1:
                raise ValueError(
                    f"{parsed.filepath.name}: expected one run_id, found {run_ids}"
                )

            df["__compile_source"] = "analysis"
            df["__compile_source_note"] = self._source_note(
                "analysis",
                parsed.filepath,
            )
            df["__compile_source_rank"] = 1
            df["__compile_source_order"] = run_order
            df["__compile_source_mtime"] = parsed.filepath.stat().st_mtime

            compile_frames.append(df)
            runinfo_by_run_id[run_ids[0]] = parsed.runinfo

        # Fail loudly if there's no data. Duh.
        if not compile_frames:
            raise ValueError(
                "No compile input data found. Provide at least one parsed workbook "
                "or a non-empty base CSV."
            )

        return compile_frames, runinfo_by_run_id


    def _source_note(self, source: str, path: str | os.PathLike[str] | None) -> str:
        """
        Return a compact provenance label for a compile input source.
        """
        if path is None or str(path).strip() == "":
            return source
        return f"{source}: {Path(path).name}"


    def _dedup_sample(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Deduplicate compile inputs by sample_id, preferring analysis over csv,
        then newer mtime, then later run order.
        """
        rank = [
            "__compile_source_rank",
            "__compile_source_mtime",
            "__compile_source_order",
        ]

        win = (
            df[["sample_id", *rank]]
            .drop_duplicates()
            .sort_values(["sample_id", *rank], ascending=True, kind="stable")
            .drop_duplicates(subset=["sample_id"], keep="last")
        )

        win_rank = df[["sample_id"]].merge(win, on="sample_id", how="left")[rank]
        keep = df[rank].eq(win_rank.to_numpy()).all(axis=1)

        if (~keep).any():
            print("- dropping duplicate inputs...")
            for sid, group in df.loc[~keep].groupby("sample_id", sort=False):
                keep_rows = df.loc[
                    keep & df["sample_id"].eq(sid),
                    ["__compile_source", "__compile_source_note"],
                ].drop_duplicates()
                drop_rows = group[
                    ["__compile_source", "__compile_source_note"]
                ].drop_duplicates()

                for _, row in keep_rows.iterrows():
                    print(
                        f"  keeping: {sid} | source={row['__compile_source']} "
                        f"| {row['__compile_source_note']}"
                    )
                for _, row in drop_rows.iterrows():
                    print(
                        f"  dropped: {sid} | source={row['__compile_source']} "
                        f"| {row['__compile_source_note']}"
                    )

        deduped = df.merge(win, on=["sample_id", *rank], how="inner")

        meta_cols = [
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
        blocks = deduped[
            [c for c in meta_cols if c in deduped.columns]
        ].drop_duplicates()
        dups = blocks.loc[blocks["sample_id"].duplicated(keep=False)]
        if not dups.empty:
            raise ValueError(
                "Duplicate sample_id values after compile dedup:\n"
                f"{dups.to_string(index=False)}"
            )

        return deduped.drop(
            columns=["__compile_source", *rank],
            errors="ignore",
        )


    def _assemble_final_combined_dataframe(self) -> pd.DataFrame:
        """
        Sort the deduplicated long dataframe into stable output order.
        """
        print("- assembling final combined dataframe...")
        sort_cols = ["samp_group", "samp_name", "samp_date"]
        return (
            self.combined_df[
                sort_cols + [c for c in self.combined_df.columns if c not in sort_cols]
            ]
            .sort_values(sort_cols, ascending=True, kind="stable")
            .reset_index(drop=True)
        )


    def _core_qc(self) -> pd.DataFrame:
        """Core QC steps applied to the combined compile dataframe."""
        df = self.combined_df.copy()

        # Re-establish numeric columns after CSV/XLSX round-trips.
        df[["bc_count", "proportion", "input"]] = (
            df[["bc_count", "proportion", "input"]]
            .apply(pd.to_numeric, errors="coerce")
            .fillna(0)
        )

        df = qc.rank_uniques(df=df)
        df["bc_name"] = df["bc_name"].str.replace(
            "SIVmac293M2", "SIVmac239M2", regex=False
        )
        df = qc.flag_shorts(df=df)

        df["ldist_samp_lvl"] = (
            df["ldist_samp_lvl"]
            .str.replace(r"^(1-)|(-1)$", "", regex=True)
            .str.replace(r"^(Unique)", "", regex=True)
            .replace("NA-NA", None)
        )

        df["multi_idx"] = df.groupby(["run_number", "samp_group", "bc_name"])[
            "samp_name"
        ].transform(lambda s: s.nunique() > 1)

        # Annotate x_idx across the whole series to catch/refresh missing or wrong in old analyses.
        df = qc.annotate_x_idx(df)

        return df


    def _build_runseries(self, source: str = "unknown") -> RunSeries:
        """
        Materialize a RunSeries from a normalized compile-mode long dataframe.
        """
        if self.combined_df.empty:
            raise ValueError(
                "Compile container construction requires a non-empty dataframe."
            )

        self.runinfo_by_run_id = {str(k): v for k, v in self.runinfo_by_run_id.items()}
        runs = {
            str(run_id): self._build_seqrun(
                run_df,
                source=source,
                runinfo=self.runinfo_by_run_id.get(str(run_id)),
            )
            for run_id, run_df in self.combined_df.groupby("run_id", sort=False)
        }

        return RunSeries(
            runs=runs,
            series_meta={"groupby_col": "run_id"},
            data_df=self.combined_df,
            runinfo_by_run_id=self.runinfo_by_run_id,
            source=source,
        )


    def _build_seqrun(
        self,
        run_df: pd.DataFrame,
        *,
        source: str,
        runinfo: RunInfo | None,
    ) -> SeqRun:
        """
        Materialize one SeqRun from a single-run slice of compiled long data.
        """
        run_ids = run_df["run_id"].dropna().astype(str).unique().tolist()
        if len(run_ids) != 1:
            raise ValueError(
                f"Compile SeqRun construction expected one run_id, found {run_ids}"
            )

        samples: dict[str, SeqSamp] = {}
        for sample_id, sample_df in run_df.groupby(
            "sample_id", dropna=False, sort=False
        ):
            sample_meta = {
                col: sample_df[col].dropna().iloc[0]
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
                ]
                if col in sample_df.columns and not sample_df[col].dropna().empty
            }
            resolved_id = str(sample_id)

            sample_meta["sample_id"] = resolved_id
            samples[resolved_id] = SeqSamp(
                df=sample_df.reset_index(drop=True),
                sample_meta=sample_meta,
                sample_id=resolved_id,
                source=source,
            )

        return SeqRun(
            samples=samples,
            run_number=self._first_value(run_df, "run_number"),
            run_name=self._first_value(run_df, "run_name"),
            run_date=self._first_value(run_df, "run_date"),
            run_id=run_ids[0],
            runinfo=runinfo,
            filepath=self._first_value(run_df, "filename"),
            source=source,
            data_df=run_df.copy(),
        )


    def _first_value(self, df: pd.DataFrame, col: str) -> Any:
        """
        Return the first non-null scalar value from a dataframe column.
        """
        if col not in df.columns or df[col].dropna().empty:
            return None
        return df[col].dropna().iloc[0]


    def _ensure_missing_runinfo(self) -> None:
        """
        Synthesize RunInfo for runs that arrived without workbook metadata.
        """
        for run_id, seq_run in self.run_series.runs.items():
            runinfo = seq_run.ensure_runinfo(
                filepath=self.runtime_settings.get("base_csv_path")
            )
            self.run_series.runinfo_by_run_id[run_id] = runinfo


    def _store_series_meta(self) -> None:
        """
        Attach compile-level metadata used by emit and interactive inspection.
        """
        self.run_series.series_meta.update(
            {
                "out_prefix": self.runtime_settings.get("out_prefix"),
                "n_source_runs": self.n_source_runs,
            }
        )

# %% .xlsx handling

class WorkbookParseError(ValueError):
    """
    Raised when a parse workbook cannot be parsed into the expected format.
    """

    def __init__(self, message: str, *, filepath: str | Path | None = None) -> None:
        """
        Store a parse failure with optional workbook context.
        """
        super().__init__(message)
        self.message = str(message)
        self.filepath = None if filepath is None else Path(filepath)

    def with_filepath(self, filepath: str | Path) -> WorkbookParseError:
        """
        Return a copy of this error bound to a specific workbook path.
        """
        return WorkbookParseError(self.message, filepath=filepath)

    def __str__(self) -> str:
        """
        Include workbook path when the error has been bound to a file.
        """
        if self.filepath is None:
            return self.message
        return f"file: {self.filepath}\nreason: {self.message}"


@dataclass
class ParsedAnalysis:
    """
    Container for one parsed Analysis or Discovery workbook.
    """

    run_id: str
    data_df: pd.DataFrame
    runinfo: RunInfo
    filepath: Path


class XlsxPathManager:
    """
    Initialize the file manager.

    Args:
        wdir: Working directory to search (defaults to the current directory).
        txtfile: Optional text file listing file paths (one per line).
        pattern: Regex pattern matching Analysis or Discovery workbooks.
    """

    def __init__(
        self,
        wdir: Optional[str] = None,
        txtfile: Optional[str] = None,
        pattern: str = r".*(?:Analysis|Discovery)\.xlsx$",
    ) -> None:
        """
        Discover workbook paths from a directory scan or explicit text file.
        """

        self.wdir: str = wdir or os.getcwd()
        self.txtfile: Optional[str] = txtfile
        self.pattern: str = pattern
        self.files: List[Path] = self._discover()

    def _discover(self) -> List[Path]:
        """
        Discover input files either from a provided txt file or by scanning the directory.
        """
        files = (
            self._fromtxt(self.txtfile) if self.txtfile else self._fromscan(self.wdir)
        )

        if not files:
            raise FileNotFoundError(
                f"No files matching pattern '{self.pattern}' found in "
                f"{'text file ' + str(self.txtfile) if self.txtfile else 'directory ' + str(self.wdir)}"
            )

        files = self._dedup_file(files)
        self.files = files
        return files

    def _fromtxt(self, txtfile: str) -> List[Path]:
        """
        Read file paths from a text file and validate their existence.
        """
        with open(txtfile, "r") as f:
            files = [Path(line.strip()) for line in f if line.strip()]

        missing = [p for p in files if not p.is_file()]
        if missing:
            mlist = "\n".join(map(str, missing))
            raise FileNotFoundError(
                f"The following input files were not found:\n{mlist}"
            )

        return files

    def _fromscan(self, directory: str) -> List[Path]:
        """
        Recursively scan the given directory for files matching the pattern.
        """
        directory_path = Path(directory)
        return [
            path
            for path in directory_path.rglob("*")
            if path.is_file()
            and not path.name.startswith("~$")
            and re.match(self.pattern, path.name)
        ]

    def _dedup_file(self, files: List[Path]) -> List[Path]:
        """
        Deduplicate by exact filename (case-sensitive), keeping the newest file.
        """
        seen: dict[str, Path] = {}
        duplicates: List[tuple[Path, Path]] = []

        for path in files:
            fname = path.name
            if fname in seen:
                newer = max([seen[fname], path], key=lambda p: p.stat().st_ctime)
                older = seen[fname] if newer == path else path
                duplicates.append((older, newer))
                seen[fname] = newer
            else:
                seen[fname] = path

        if duplicates:
            logger.warning("Duplicate filenames found; keeping the newest file:")
            for old, new in duplicates:
                logger.warning(f"  Older: {old}, Newer: {new}")

        return list(seen.values())


class AnalysisWorkbookParser:
    """
    AnalysisWorkbookParser class:
    Parses an Excel analysis file, extracts metadata and sample tables,
    and normalizes column names.
    """

    def __init__(self, filepath: str | Path) -> None:
        """
        Load and parse one Analysis or Discovery workbook into normalized data.
        """

        self.filepath: Path = Path(filepath)
        try:
            self._load_xlsx()
            self._parse_xlsx()
        except WorkbookParseError as e:
            raise e.with_filepath(self.filepath) from e
        except Exception as e:
            raise WorkbookParseError(
                f"Unexpected {type(e).__name__}: {e}",
                filepath=self.filepath,
            ) from e

    def _load_xlsx(self, **kwargs: Any) -> None:
        """
        Load all sheets from the Excel file into self.sheets.
        """
        # logger.debug(f"Loading Excel file: {self.filepath}")
        self.sheets: Dict[str, pd.DataFrame] = pd.read_excel(
            self.filepath,
            sheet_name=None,  # load all sheets into a dict
            header=None,
            engine="openpyxl",
            **kwargs,
        )

    def _parse_xlsx(self) -> None:
        """
        Process the loaded Excel workbook into a standardized long-format DataFrame.
        """
        all_samples: List[pd.DataFrame] = []

        # Obtain runinfo from first sheet
        first_sheet = next(iter(self.sheets.values()))
        self.runinfo = RunInfo.from_table(
            raw_df=self._sheet_to_runinfo_df(first_sheet),
            filepath=str(self.filepath),
        )
        # logger.debug(f"Runinfo: {self.runinfo.run_id}")

        # Process only 'samples' sheets (skip first and non-sample sheets)
        for i, (sname, df) in enumerate(self.sheets.items()):
            if i == 0:
                continue
            if "samples" not in sname.lower():
                continue

            # logger.debug(f"Parsing: {sname}")
            samples = self._parse_samples(df)
            all_samples.extend(samples)

        self.run_id = make_run_id(
            self.runinfo.run_number,
            self.filepath.name,
            self.filepath,
        )

        self.data_df = self._combine_sample_tables(all_samples)
        self.analysis = ParsedAnalysis(
            run_id=self.run_id,
            data_df=self.data_df.copy(),
            runinfo=self.runinfo,
            filepath=self.filepath,
        )

    def _sheet_to_runinfo_df(self, first_df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert the first sheet of a parse workbook back into a headered runinfo table.
        """
        if first_df.empty:
            return pd.DataFrame()

        header = first_df.iloc[0].fillna("").astype(str).str.strip()
        runinfo_df = first_df.iloc[1:].reset_index(drop=True).copy()
        runinfo_df.columns = header
        runinfo_df = runinfo_df.loc[:, runinfo_df.columns != ""]
        return runinfo_df

    def _parse_samples(self, df: pd.DataFrame) -> List[pd.DataFrame]:
        """
        Extract multiple side-by-side tables from a 'samples' sheet DataFrame.
        """
        r, c = self._detect_breaks(df)

        # Recognizes specific sheet formats:
        valid_breaks = {
            (3, 6),  # r=3, c=6: original missing 'short bcs'
            (3, 7),  # r=3, c=7: original 7 data columns
            (7, 7),  # r=7, c=7: current non-legacy analysis output
            (8, 7),  # r=8, c=7: pre-refactor non-legacy analysis output
        }
        if r is None or c is None or (r, c) not in valid_breaks:
            raise WorkbookParseError(
                f"Nonstandard format in samples sheet (r={r}, c={c})."
            )

        header_row = df.iloc[r + 1]

        # Define data column blocks
        mask = header_row.notna()
        block_ids = (mask != mask.shift(fill_value=False)).cumsum() * mask
        col_blocks = [
            group.index.tolist() for _, group in header_row[mask].groupby(block_ids)
        ]

        tables: List[pd.DataFrame] = []

        for block in col_blocks:
            meta_df = df.iloc[0:r, block]
            meta_dict = self._extract_metadata(meta_df, r)

            block_df = df.iloc[r + 2 :, block].copy()
            block_df.columns = df.iloc[r + 1, block]

            # Define above_cutoff and drop spacer rows
            is_blank = block_df.isna().all(axis=1)
            block_df["above_cutoff"] = (
                is_blank.cumsum() == 0
            )  # this pulls first block, non-derivative and > 1/input
            block_df = block_df[~is_blank].reset_index(drop=True)

            sample_meta = self._resolve_runinfo_sample(meta_dict)
            for k, v in sample_meta.items():
                block_df[k] = v

            block_df["run_number"] = self.runinfo.run_number

            # Add shorts col if missing
            if c == 6:
                block_df["Short barcode?"] = "not checked!"

            # Standardize column names per sample block so each table remains one SeqSamp.
            block_df = block_df.rename(
                columns={
                    "Barcode": "bc_name",
                    "Sequence": "bc_seq",
                    "Counts": "bc_count",
                    "Proportion": "proportion",
                    "Hamming Dist to Another Sequence": "ldist_samp_lvl",
                    "putative_parent": "ldist_samp_lvl",
                    "Index hopping?": "x_idx",
                    "Multi-group share": "x_idx",
                    "Multi-group?": "x_idx",
                    "x_group": "x_idx",
                    "Short barcode?": "short_bc",
                }
            )

            if "x_idx" in block_df.columns:
                block_df["x_idx"] = pd.to_numeric(block_df["x_idx"], errors="coerce")

            tables.append(block_df.reset_index(drop=True))

        return tables

    def _detect_breaks(self, df: pd.DataFrame) -> Tuple[Optional[int], Optional[int]]:
        """
        Detect header break row (r) and determine number of data columns (c).

        Returns:
            Tuple of (r, c) where:
            - r is the row index of the break
            - c is the count of header columns
        """
        mask = df.isna().all(axis=1).to_numpy()
        r = int(mask.argmax()) if mask.any() else None

        if r is not None:
            header_row = df.iloc[r + 1]
            if header_row.isna().any():
                c = int(header_row.isna().to_numpy().argmax())
            else:
                c = int(header_row.notna().sum())
        else:
            c = None

        return r, c

    def _extract_metadata(self, meta_block: pd.DataFrame, r: int) -> Dict[str, Any]:
        """
        Extract sample metadata from a block given the detected break row.
        """
        if r == 3:
            return {
                "samp_group": meta_block.iat[0, 0],
                "samp_name": meta_block.iat[0, 1],
                "samp_date": normalize_date(meta_block.iat[0, 2]),
                "input": meta_block.iat[2, 0],
            }
        if r == 7:
            return {
                "samp_group": meta_block.iat[3, 1],
                "samp_name": meta_block.iat[4, 1],
                "samp_date": normalize_date(meta_block.iat[5, 1]),
                "input": meta_block.iat[6, 1],
            }
        if r == 8:
            return {
                "samp_group": meta_block.iat[4, 1],
                "samp_name": meta_block.iat[5, 1],
                "samp_date": normalize_date(meta_block.iat[6, 1]),
                "input": meta_block.iat[7, 1],
            }
        raise ValueError(f"Unrecognized header break row r={r}")

    def _resolve_runinfo_sample(
        self,
        meta_dict: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Recover runinfo-derived sample metadata for a sample block.
        """
        match_df = self.runinfo.sample_rows
        for key in ["samp_group", "samp_name", "samp_date"]:
            match_df = match_df.loc[match_df[key] == meta_dict.get(key)]

        if len(match_df) != 1:
            raise WorkbookParseError(
                "Could not uniquely match sample sheet metadata to runinfo."
            )

        row = match_df.iloc[0]
        return {str(k): v for k, v in row.to_dict().items()}

    def _combine_sample_tables(self, sample_tables: List[pd.DataFrame]) -> pd.DataFrame:
        """
        Concatenate parsed sample blocks into one normalized long-format dataframe.
        """
        if sample_tables:
            df = pd.concat([table.copy() for table in sample_tables], ignore_index=True)
        else:
            df = pd.DataFrame()

        df["filename"] = self.filepath.name
        df["run_id"] = self.run_id

        df = normalize_long_df(
            df,
            run_number=self.runinfo.run_number,
            run_name=self.runinfo.run_name,
            run_date=self.runinfo.run_date,
            source_path=self.filepath,
        )

        df["run_id"] = df["run_id"].where(df["run_id"].notna(), self.run_id)

        return df

# %% Versions
"""
v1.0.0 20260625
    - Created as the linear compile-mode pipeline module.
    - Combines former ingest/compile.py source normalization, sample-level
      deduplication, compile QC, RunSeries construction, and missing RunInfo
      synthesis machinery.
    - Combines former ingest/xlsx.py parse-workbook discovery and parsing
      machinery.
    - Keeps workbook/csv emission outside CompileMode in __main__/emit.
"""
