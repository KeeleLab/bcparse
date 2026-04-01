# runinfo.py
"""
Name:       runinfo.py
Author:     CAG
Version:    1.0
Date:       2026/03/26

Normalized metadata container for one weirdly shaped but important runinfo table.
"""

#%% Imports

from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

#%% Global helpers

def make_sample_id(
    samp_group: Any,
    samp_name: Any,
    samp_date: Any,
    run_number: Any,
) -> str:
    """
    Generate the composite sample identity used by compile mode.
    """
    date_val = pd.to_datetime(samp_date, errors="coerce")
    date_str = "" if pd.isna(date_val) else date_val.strftime("%Y-%m-%d")

    return "::".join(
        [
            "" if pd.isna(samp_group) else str(samp_group),
            "" if pd.isna(samp_name) else str(samp_name),
            date_str,
            "" if pd.isna(run_number) else str(run_number),
        ]
    )

def make_run_id(
    run_number: Any,
    filename: Any = None,
    run_name: Any = None,
) -> str | None:
    """
    Generate a stable run_id used to tie long-format data back to a RunInfo.
    """
    run_number_str = "" if pd.isna(run_number) else str(run_number)
    filename_str = "" if pd.isna(filename) else str(filename)
    run_name_str = "" if pd.isna(run_name) else str(run_name)

    if filename_str:
        return "::".join([run_number_str, filename_str])

    if run_name_str:
        return "::".join([run_number_str, run_name_str])

    if run_number_str:
        return run_number_str

    return None


#%% RunInfo model

@dataclass
class RunInfo:
    """
    Normalized metadata companion for a single sequencing run.

    This object carries both sample identifiers:
    - `idx_name`: parse-mode sample identity derived from runinfo index columns
    - `sample_id`: compile-mode composite identity derived from biological metadata
    """
    run_meta: dict[str, Any] = field(default_factory=dict)                # Run-level metadata, e.g. run_number, run_name, run_date, filename.
    sample_meta: dict[str, dict[str, Any]] = field(default_factory=dict)  # Sample metadata keyed by sample identity; in normal parse/workbook flows that key is idx_name.
    raw_df: pd.DataFrame | None = None
    source: str = "unknown"
    filepath: Path | None = None

    def __post_init__(self) -> None:
        # Make private copies of mutable inputs.
        self.run_meta = dict(self.run_meta)
        self.sample_meta = {
            str(sample_id): dict(meta)
            for sample_id, meta in self.sample_meta.items()
        }
        self.raw_df = None if self.raw_df is None else self.raw_df.copy()

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

        # Fill the compile-mode composite sample id for each sample.
        run_number = self.run_meta.get("run_number")
        for meta in self.sample_meta.values():
            meta.setdefault(
                "sample_id",
                make_sample_id(
                    meta.get("samp_group"),
                    meta.get("samp_name"),
                    meta.get("samp_date"),
                    run_number,
                ),
            )

    # ====================
    # INGEST
    # ====================

    @classmethod
    def from_table(
        cls,
        *,
        raw_runinfo: pd.DataFrame,
        filepath: str | None = None,
        dual_index: bool | None = None,
    ) -> "RunInfo":
        """
        Normalize a raw runinfo workbook table into a RunInfo container.
        """
        # Clean the incoming workbook table first.
        runinfo = raw_runinfo.copy()
        runinfo.columns = runinfo.columns.str.strip()
        runinfo = runinfo.apply(
            lambda col: col.apply(lambda x: x.strip() if isinstance(x, str) else x)
        )

        if "Date" in runinfo.columns:
            runinfo["Date"] = pd.to_datetime(runinfo["Date"], errors="coerce").dt.normalize()

        if "Input TOTAL PER BARCODE" in runinfo.columns:
            runinfo["Input TOTAL PER BARCODE"] = cls._normalize_input_col(
                runinfo["Input TOTAL PER BARCODE"]
            )

        # Pull run-level metadata from the legacy runinfo layout.
        run_number = cls._extract_run_number(runinfo)
        run_meta = {
            "run_number": run_number,
            "run_name": cls._label_value(runinfo, "Run Name"),
            "run_date": cls._label_value(runinfo, "Date"),
        }

        if filepath is not None:
            run_meta["filepath"] = filepath
            run_meta["filename"] = os.path.basename(filepath)

        # Build canonical sample metadata from sample-bearing rows.
        sample_meta = cls._build_sample_meta(runinfo, run_number, dual_index=dual_index)

        return cls(
            run_meta=run_meta,
            sample_meta=sample_meta,
            raw_df=runinfo,
            source="runinfo",
            filepath=filepath,
        )

    # ====================
    # COMPILE
    # ====================

    @staticmethod
    def _normalize_input_col(series: pd.Series) -> pd.Series:
        """
        Validate and coerce the runinfo input column to numeric/inf values.
        """
        series = series.replace(
            to_replace=r"(?i)^inf$",
            value=np.inf,
            regex=True,
        )

        invalid = (
            ~series.astype(str)
            .str.strip()
            .str.match(r"^(inf|nan|\d+(\.\d+)?|\.\d+)$", na=False)
        )

        if invalid.any():
            invalid_strs = [str(x) for x in series.loc[invalid].unique()]
            raise ValueError(f"Invalid input: {invalid_strs} - use numeric or inf.")

        return pd.to_numeric(series, errors="raise")


    @staticmethod
    def _extract_run_number(runinfo: pd.DataFrame) -> Any:
        """
        Pull the run number from the legacy runinfo layout.
        """
        if "Run Number" not in runinfo.columns:
            return None

        run_number_col = runinfo["Run Number"]
        return run_number_col.iloc[0] if not run_number_col.empty else None


    @staticmethod
    def _label_value(runinfo: pd.DataFrame, label: str) -> Any:
        """
        Read add'l labeled value from the legacy runinfo 1st column.
        """
        if "Run Number" not in runinfo.columns:
            return None

        matches = runinfo.index[runinfo["Run Number"] == label].tolist()
        if not matches:
            return None

        next_idx = matches[0] + 1
        if next_idx >= len(runinfo.index):
            return None

        return runinfo.loc[next_idx, "Run Number"]


    @classmethod
    def _build_sample_meta(
        cls,
        runinfo: pd.DataFrame,
        run_number: Any,
        *,
        dual_index: bool | None = None,
    ) -> dict[str, dict[str, Any]]:
        """
        Build per-sample metadata keyed by idx_name.
        """
        sample_meta: dict[str, dict[str, Any]] = {}

        required_cols = {"Barcodes", "Date", "Input TOTAL PER BARCODE"}
        if not required_cols.issubset(runinfo.columns):
            return sample_meta

        sample_rows = runinfo.dropna(subset=["Barcodes"]).copy()
        if "idx_name" in sample_rows.columns and sample_rows["idx_name"].notna().any():
            sample_rows["idx_name"] = sample_rows["idx_name"].astype(str).str.strip()
            runinfo.loc[sample_rows.index, "idx_name"] = sample_rows["idx_name"]
        else:
            sample_rows["idx_name"] = sample_rows["Barcodes"]
            runinfo["idx_name"] = runinfo.get("Barcodes")

            if dual_index and "(F Barcode)" in sample_rows.columns:
                has_p7 = sample_rows["(F Barcode)"].notna() & sample_rows["(F Barcode)"].astype(str).ne("")
                sample_rows.loc[has_p7, "idx_name"] = (
                    sample_rows.loc[has_p7, "Barcodes"].astype(str)
                    + "_"
                    + sample_rows.loc[has_p7, "(F Barcode)"].astype(str)
                )

            runinfo.loc[sample_rows.index, "idx_name"] = sample_rows["idx_name"]

        keep_cols = [
            colname
            for colname in ["Animal", "Sample", "idx_name", "Date", "Input TOTAL PER BARCODE"]
            if colname in sample_rows.columns
        ]

        sample_rows = sample_rows[keep_cols].rename(
            columns={
                "Animal": "samp_group",
                "Sample": "samp_name",
                "Date": "samp_date",
                "Input TOTAL PER BARCODE": "input",
            }
        )

        for _, row in sample_rows.iterrows():
            idx_name = row.get("idx_name")
            if pd.isna(idx_name):
                continue

            meta = {k: v for k, v in row.to_dict().items() if pd.notna(v)}
            meta["sample_id"] = make_sample_id(
                row.get("samp_group"),
                row.get("samp_name"),
                row.get("samp_date"),
                run_number,
            )
            sample_meta[str(row.get("idx_name"))] = meta

        return sample_meta

    # ====================
    # HELPERS
    # ====================

    def find_idx_sample_id_collisions(self) -> dict[str, list[str]]:
        """
        Find cases where distinct runinfo index labels collapse to the same
        compile-mode composite sample identity.
        """
        sample_meta_df = self.sample_meta_df
        if sample_meta_df.empty:
            return {}

        required_cols = {"idx_name", "sample_id"}
        if not required_cols.issubset(sample_meta_df.columns):
            return {}

        sample_rows = sample_meta_df.dropna(subset=["idx_name", "sample_id"]).copy()
        if sample_rows.empty:
            return {}

        grouped = (
            sample_rows.groupby("sample_id")["idx_name"]
            .apply(lambda s: sorted({str(value) for value in s.dropna().tolist() if str(value)}))
            .to_dict()
        )

        return {
            sample_id: idx_names
            for sample_id, idx_names in grouped.items()
            if len(idx_names) > 1
        }

    def get_sample_meta_by_idx_name(self, idx_name: str, default: Any = None) -> Any:
        for meta in self.sample_meta.values():
            if str(meta.get("idx_name")) == str(idx_name):
                return meta
        return default

    # ====================
    # EMIT
    # ====================

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
    def sample_meta_df(self) -> pd.DataFrame:
        if not self.sample_meta:
            return pd.DataFrame()
        return pd.DataFrame(self.sample_meta.values())

    @property
    def full_index_col(self) -> str:
        """
        Return the canonical index label column available in raw runinfo.
        """
        if self.raw_df is None:
            return "idx_name"
        if "idx_name" in self.raw_df.columns:
            return "idx_name"
        if "Barcodes" in self.raw_df.columns:
            return "Barcodes"
        return "idx_name"

    @property
    def expected_groups(self) -> list[str]:
        """
        Return the expected sample groups declared by metadata/runinfo.
        """
        if self.raw_df is not None and "Animal" in self.raw_df.columns:
            return self.raw_df["Animal"].dropna().astype(str).unique().tolist()

        sample_meta_df = self.sample_meta_df
        if sample_meta_df.empty or "samp_group" not in sample_meta_df.columns:
            return []

        return sample_meta_df["samp_group"].dropna().astype(str).unique().tolist()

    def copy_raw_df(self) -> pd.DataFrame:
        return pd.DataFrame() if self.raw_df is None else self.raw_df.copy()

#%% Versions
"""
v1.0
- Initial version, accomodating new model/container structure
- New dedicated class to manage runinfo-associated metadata
"""
