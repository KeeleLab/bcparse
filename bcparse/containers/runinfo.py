# runinfo.py
"""
Name:      runinfo.py
Author:     CAG
Version:    2.0.1
Date:       2026/06/23

Metadata container for one sequencing run's runinfo Excel table.

Layout contract (columns):
    Run Number | Animal | Sample | Date | Barcodes | Input TOTAL PER BARCODE |
    Input TOTAL PER WELL | (F Barcode) | (cDNA) | [passthrough cols...]

Run-level metadata lives in the first column below the sample block:
    "Run Name"  → value in next row
    "Date"      → value in next row
    Free-text notes below these are ignored.
"""

# %% Imports

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from bcparse.qc.identity import (
    make_run_id,
    make_sample_id,
    normalize_date,
    normalize_label,
)

# %% Helpers

# Expected cols from input runinfo.
_INTERNAL_COLS = {
    "Run Number",  # column 0, row 0, contains other metadata labels in subsequent rows
    "Animal",  # Req'd -> samp_group
    "Sample",  # Req'd -> samp_name
    "Date",  # Req'd -> samp_date
    "Barcodes",  # Req'd -> idx_name
    "Input TOTAL PER BARCODE",  # Req'd -> input
    "Input TOTAL PER WELL",  # legacy, optional
    "(F Barcode)",  # legacy, optional
    "(cDNA)",  # legacy, optional
}

# %% RunInfo


class RunInfo:
    """
    Parsed metadata for one sequencing run.

    Attributes
    ----------
    run_number      : value from row 1 of the Run Number column
    run_name        : value from the "Run Name" label block
    run_date        : value from the "Date" label block
    filename        : source file name
    run_id          : stable composite run identifier
    sample_rows     : normalized per-sample DataFrame (one row per index)
    raw_df          : the original table as read, for output round-trips
    passthrough_cols: column names carried forward into long-format data
    """

    def __init__(
        self,
        run_number,
        run_name,
        run_date,
        sample_rows: pd.DataFrame,
        raw_df: pd.DataFrame,
        passthrough_cols: list[str],
        filename: str | None = None,
        filepath: str | Path | None = None,  # parse -> path, compile -> str | None
        dual_index: bool = False,
    ):
        self.run_number = run_number
        self.run_name = run_name
        self.run_date = normalize_date(run_date)
        self.filename = filename
        self.filepath = Path(filepath) if filepath else None
        self.run_id = make_run_id(run_number, filename, filepath)
        self.sample_rows = sample_rows  # normalized; one row per index
        self.raw_df = raw_df.copy()
        self.passthrough_cols = passthrough_cols
        self.dual_index = dual_index

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_table(
        cls,
        raw_df: pd.DataFrame,
        filepath: str | Path | None = None,
        dual_index: bool = False,
    ) -> "RunInfo":
        """
        Parse a raw runinfo Excel table into a RunInfo object.
        Called from:
            - ingest/parse_mode.py
            - ingest/compile_mode.py, for Analysis/Discovery workbook first sheets
        """

        # --- Preprocess, strip whitespace ---
        df = raw_df.copy()
        df.columns = df.columns.str.strip()
        df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

        # --- Run-level metadata ---
        # Run Number is in row 0 of the "Run Number" column.
        run_number = df["Run Number"].iloc[0] if "Run Number" in df.columns else None
        run_number = None if pd.isna(run_number) else run_number

        # Run Name and Date live as label→next-row pairs in the Run Number column.
        run_name = cls._label_value(df, "Run Name")
        run_date = cls._label_value(df, "Date")

        # --- Sample rows ---
        # Any row with a non-null `Barcodes` cell is a sample row.
        sample_mask = (
            df["Barcodes"].notna()
            if "Barcodes" in df.columns
            else pd.Series(False, index=df.index)
        )
        sample_df = df.loc[sample_mask].copy()

        # Build idx_name per row
        # P5 only,
        sample_df["idx_name"] = sample_df["Barcodes"].astype(str)
        # or "P5_P7" for dual-index runs.
        if dual_index and "(F Barcode)" in sample_df.columns:
            has_p7 = sample_df["(F Barcode)"].notna() & sample_df["(F Barcode)"].astype(
                str
            ).str.strip().ne("")
            sample_df.loc[has_p7, "idx_name"] = (
                sample_df.loc[has_p7, "Barcodes"].astype(str)
                + "_"
                + sample_df.loc[has_p7, "(F Barcode)"].astype(str)
            )

        # Rename columns to internal names
        sample_df = sample_df.rename(
            columns={
                "Animal": "samp_group",
                "Sample": "samp_name",
                "Date": "samp_date",
                "Input TOTAL PER BARCODE": "input",
                "(F Barcode)": "p7_name",
            }
        )

        # Normalize values
        sample_df["samp_date"] = sample_df["samp_date"].apply(normalize_date)
        sample_df["input"] = cls._normalize_input(sample_df["input"])
        sample_df["samp_group"] = (
            sample_df["samp_group"].map(normalize_label).astype("string")
        )

        # Generate a stable sample_id for each row
        sample_df["sample_id"] = sample_df.apply(
            lambda r: make_sample_id(
                r["samp_group"], r["samp_name"], r["samp_date"], run_number, filepath
            ),
            axis=1,
        )

        # Passthrough: any column not in the internal set.
        # Allows user-defined columns to be carried forward into long-format data.
        passthrough_cols = [
            c
            for c in sample_df.columns
            if c not in _INTERNAL_COLS
            and c
            not in {
                "samp_group",
                "samp_name",
                "samp_date",
                "input",
                "p7_name",
                "idx_name",
                "sample_id",
                "Barcodes",
            }
        ]

        # Keep req'd internal columns and user-supplied passthrough columns for sample_rows.
        keep_cols = [
            "idx_name",
            "samp_group",
            "samp_name",
            "samp_date",
            "input",
            "p7_name",
            "sample_id",
        ] + passthrough_cols

        # Build the final sample_rows DataFrame with only the kept columns, in order.
        sample_rows = sample_df[
            [c for c in keep_cols if c in sample_df.columns]
        ].reset_index(drop=True)

        # Pull filename from filepath
        fn = Path(filepath).name if filepath else None

        # Return the constructed RunInfo object
        return cls(
            run_number=run_number,
            run_name=run_name,
            run_date=run_date,
            sample_rows=sample_rows,
            raw_df=raw_df,
            passthrough_cols=passthrough_cols,
            filename=fn,
            filepath=filepath,
            dual_index=dual_index,
        )

    @classmethod
    def from_long_df(
        cls,
        run_df: pd.DataFrame,
        *,
        filepath: str | Path | None = None,
    ) -> "RunInfo":
        """
        Synthesize a minimal RunInfo from one run-scoped long-format dataframe.

        Called from:
            - runseries.py slices compiled data to seqruns
            - This is called via SeqRun.ensure_runinfo()
        """

        # --- Validate required columns ---
        df = run_df.copy()

        required_cols = [
            "run_number",
            "run_name",
            "run_date",
            "filename",
            "sample_id",
            "idx_name",
            "samp_group",
            "samp_name",
            "samp_date",
            "input",
        ]

        missing = [col for col in required_cols if col not in df.columns]

        if missing:
            raise ValueError(
                f"RunInfo.from_long_df: missing required column(s): {missing}"
            )

        # --- Normalize values ---
        df["run_date"] = df["run_date"].apply(normalize_date)
        df["samp_date"] = df["samp_date"].apply(normalize_date)
        df["input"] = pd.to_numeric(df["input"], errors="coerce")

        # --- Extract run metadata ---
        run_number = cls._single_value(df["run_number"], "run_number")
        run_name = cls._single_value(df["run_name"], "run_name")
        run_date = cls._single_value(df["run_date"], "run_date")
        filename = cls._single_value(df["filename"], "filename")

        if filename is None and filepath is not None:
            filename = Path(filepath).name

        # --- Get discrete sample records per run ---
        sample_records = []

        for sample_id, sample_df in df.groupby("sample_id", dropna=False, sort=False):
            idx_name = cls._single_value(
                sample_df["idx_name"], "idx_name", sample_id=sample_id
            )
            samp_group = cls._single_value(
                sample_df["samp_group"], "samp_group", sample_id=sample_id
            )
            samp_name = cls._single_value(
                sample_df["samp_name"], "samp_name", sample_id=sample_id
            )
            samp_date = cls._single_value(
                sample_df["samp_date"], "samp_date", sample_id=sample_id
            )
            input_val = cls._single_value(
                sample_df["input"], "input", sample_id=sample_id
            )
            barcode, p7_name = cls._split_idx_name(idx_name)

            sample_records.append(
                {
                    "Run Number": run_number,
                    "Animal": samp_group,
                    "Sample": samp_name,
                    "Date": samp_date,
                    "Barcodes": barcode,
                    "Input TOTAL PER BARCODE": input_val,
                    "Input TOTAL PER WELL": "",
                    "(F Barcode)": "" if p7_name is None else p7_name,
                    "(cDNA)": "",
                    "full_idx": idx_name,
                    "notes": pd.NA,
                    "sample_id": sample_id,
                }
            )

        # --- Build raw_df with legacy names ---
        raw_df = pd.DataFrame(
            sample_records,
            columns=[
                "Run Number",
                "Animal",
                "Sample",
                "Date",
                "Barcodes",
                "Input TOTAL PER BARCODE",
                "Input TOTAL PER WELL",
                "(F Barcode)",
                "(cDNA)",
                "full_idx",
                "notes",
                "sample_id",
            ],
        )

        # --- Build sample_rows with normalized names ---
        sample_rows = raw_df.rename(
            columns={
                "Animal": "samp_group",
                "Sample": "samp_name",
                "Date": "samp_date",
                "Input TOTAL PER BARCODE": "input",
                "(F Barcode)": "p7_name",
                "full_idx": "idx_name",
            }
        )

        sample_rows["p7_name"] = sample_rows["p7_name"].replace("", None)

        sample_rows = sample_rows[
            [
                "idx_name",
                "samp_group",
                "samp_name",
                "samp_date",
                "input",
                "p7_name",
                "notes",
                "sample_id",
            ]
        ].reset_index(drop=True)

        dual_index = bool(
            raw_df["full_idx"].astype(str).str.contains("_", na=False).any()
        )

        # --- Return the constructed RunInfo object ---
        return cls(
            run_number=run_number,
            run_name=run_name,
            run_date=run_date,
            sample_rows=sample_rows,
            raw_df=raw_df,
            passthrough_cols=["full_idx", "notes", "sample_id"],
            filename=filename,
            filepath=filepath,
            dual_index=dual_index,
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _single_value(series: pd.Series, field: str, sample_id=None):
        """
        Return the single unique value in a series, or raise ValueError if multiple.
        """
        vals = series.dropna().drop_duplicates().tolist()
        if len(vals) > 1:
            if sample_id is None:
                raise ValueError(
                    f"RunInfo.from_long_df: '{field}' has conflicting values within run: {vals}"
                )
            raise ValueError(
                f"RunInfo.from_long_df: '{field}' has conflicting values for sample_id={sample_id}: {vals}"
            )
        return None if not vals else vals[0]

    @staticmethod
    def _split_idx_name(idx_name):
        """
        Split an idx_name into P5 and P7 components if present.
        """
        if pd.isna(idx_name):
            return None, None
        parts = str(idx_name).split("_", 1)
        if len(parts) == 1:
            return parts[0], None
        return parts[0], parts[1]

    @staticmethod
    def _label_value(df: pd.DataFrame, label: str):
        """
        Return the cell immediately below the first occurrence of `label` in column 0.
        """
        col = df.iloc[:, 0]
        matches = col.index[col.astype(str).str.strip() == label].tolist()
        if not matches:
            return None
        next_idx = matches[0] + 1
        return df.iloc[:, 0].iloc[next_idx] if next_idx < len(df) else None

    @staticmethod
    def _normalize_input(series: pd.Series) -> pd.Series:
        """
        Coerce input column to numeric, accepting 'inf' as a valid value.
        """
        s = series.astype(str).str.strip().str.replace(",", "", regex=False)
        s = s.str.replace(r"(?i)^inf$", "inf", regex=True)
        invalid = ~s.str.match(r"^(inf|nan|\d+(\.\d+)?|\.\d+)$", na=False)
        if invalid.any():
            raise ValueError(
                f"Invalid input values: {series[invalid].unique().tolist()}"
            )
        return pd.to_numeric(s, errors="coerce").replace(float("inf"), np.inf)

    # ------------------------------------------------------------------
    # Accessors used downstream
    # ------------------------------------------------------------------

    def index_key_for_output(self, runinfo_df: pd.DataFrame) -> pd.Series:
        """
        Return the barcode index key for a runinfo-shaped output table.
        """
        key = runinfo_df["Barcodes"].copy()
        if self.dual_index:
            has_p7 = (
                key.notna()
                & runinfo_df["(F Barcode)"].notna()
                & runinfo_df["(F Barcode)"].astype(str).str.strip().ne("")
            )
            key.loc[has_p7] = (
                key.loc[has_p7].astype(str)
                + "_"
                + runinfo_df.loc[has_p7, "(F Barcode)"].astype(str)
            )
        return key

    @property
    def expected_groups(self) -> list[str]:
        """
        Return the distinct sample groups in this runinfo, or an empty list if none.
        Used in emit/workbooks.py to write the "Expected Groups" sheet in parse workbooks.
        """
        if "samp_group" not in self.sample_rows.columns:
            return []
        return self.sample_rows["samp_group"].dropna().astype(str).unique().tolist()

    def find_idx_sample_id_collisions(self) -> dict[str, list[str]]:
        """
        Find distinct idx_name values that map to the same composite sample_id.
        Used in main.py before parsing to warn users of potential collisions.
        """
        if self.sample_rows.empty or not {"idx_name", "sample_id"}.issubset(
            self.sample_rows.columns
        ):
            return {}
        grouped = (
            self.sample_rows.groupby("sample_id")["idx_name"]
            .apply(lambda s: sorted(s.dropna().astype(str).unique().tolist()))
            .to_dict()
        )
        return {str(sid): list(idxs) for sid, idxs in grouped.items() if len(idxs) > 1}

    def copy_raw_df(self) -> pd.DataFrame:
        """
        Return a copy of the raw DataFrame.
        Used in emit/workbooks.py to write runinfo tables to output workbooks.
        """
        return self.raw_df.copy()


# %%  Repository usage notes
"""
RunInfo is the run-level metadata container used throughout parse and
compile flows.

Construction:
 - Parse mode reads runinfo.xlsx in ingest/parse_mode.py via RunInfo.from_table().
 - Compile mode reconstructs RunInfo from parse-workbook first sheets in
   ingest/compile_mode.py via RunInfo.from_table().
 - SeqRun.ensure_runinfo() synthesizes RunInfo from long-format data via
   RunInfo.from_long_df() in compiled/base CSV inputs.

Storage and flow:
 - ParseMode and ParsedAnalysis both carry RunInfo during ingest.
 - SeqRun stores one optional RunInfo as seq_run.runinfo.
 - RunSeries stores per-run metadata in runinfo_by_run_id.

Downstream consumers:
 - main checks find_idx_sample_id_collisions() before parsing FASTQ.
 - ParseMode uses raw_df to validate expected P5 or P5/P7 indexes and
   sample_rows to merge normalized sample metadata into count data.
 - emit/workbooks.py uses copy_raw_df(), expected_groups, and
   index_key_for_output() when writing parse workbooks.
 - emit/views.py concatenates per-run raw runinfo tables for compile output.

Entry points:
 - Parse: runinfo_path is a required parse setting in CLI, GUI, and ParseSettings.
 - Compile: Pulled from parse-workbook first sheets, or synthesized from base CSVs.
"""

# %% Versions
"""
v2.0.1 20260623
 - Format and comment updates
 - Minor edits to soothe pylint warnings

v2.0.0 20260526
 - Refactored into bcparse package structure
 - Added from_long_df synthetic runinfo constructor
"""
