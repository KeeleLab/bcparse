# identity.py
"""
Name:       identity.py
Author:     CAG
Version:    3.0
Date:       2026/03/20
Refactored: 2026/05/26

Consolidating qc functions for processing barcode data.
These are used across the parse- and compile-side build pipelines.
"""

#%% Imports
from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd

#%% General

def make_run_id(run_number, filename, source_path=None) -> str:
    required = {
        "Run Number": run_number,
        "filename": filename,
    }
    missing = [k for k, v in required.items() if pd.isna(v) or str(v).strip() == ""]
    if missing:
        src = Path(source_path).name if source_path else "input"
        raise ValueError(f"{src}: missing required run identity field(s): {', '.join(missing)}")

    return "::".join([str(run_number), str(filename)])


def make_sample_id(samp_group, samp_name, samp_date, run_number, source_path=None) -> str:
    required = {
        "Animal": samp_group,
        "Sample": samp_name,
        "Date": samp_date,
        "Run Number": run_number,
    }
    missing = [k for k, v in required.items() if pd.isna(v) or str(v).strip() == ""]
    if missing:
        src = Path(source_path).name if source_path else "input"
        raise ValueError(f"{src}: missing required sample identity field(s): {', '.join(missing)}")

    return "::".join([str(samp_group), str(samp_name), str(samp_date), str(run_number)])


def normalize_date(date_val) -> Any:
    """
    Coerce a value to canonical YYYY-MM-DD format if possible.
    Preserve unparseable non-missing text.
    """
    dt = pd.to_datetime(date_val, errors="coerce")
    if pd.isna(dt):
        return date_val
    return dt.strftime("%Y-%m-%d")


def normalize_label(value) -> Any:
    """
    Coerce metadata labels to text, dropping Excel-derived .0 suffixes.
    """
    if pd.isna(value):
        return pd.NA
    if isinstance(value, (float, np.floating)):
        return str(int(value))
    return str(value).strip()


def coerce_boolish(value):
    """
    Coerce common workbook/CSV representations to bool where possible.
    """
    if pd.isna(value):
        return False

    if isinstance(value, bool):
        return value

    if isinstance(value, (int, float)):
        return bool(value)

    val = str(value).strip().lower()

    if val in {"true", "yes", "y", "1"}:
        return True
    if val in {"false", "no", "n", "0", ""}:
        return False

    return False


def normalize_long_df(
    df: pd.DataFrame,
    *,
    run_number=None,
    run_name=None,
    run_date=None,
    source_path=None,
    ) -> pd.DataFrame:
    """
    Normalize a long-format sequencing dataframe to the shared column contract.
    """
    df = df.copy()

    defaults = {
        "run_id": None,
        "run_number": run_number,
        "run_name": run_name,
        "run_date": run_date,
        "sample_id": None,
        "idx_name": None,
        "samp_group": None,
        "samp_name": None,
        "samp_date": None,
        "input": None,
        "bc_name": None,
        "bc_count": 0,
        "proportion": 0.0,
        "above_cutoff": False,
        "x_idx": np.nan,
        "multi_idx": False,
        "short_bc": False,
        "putative_parent": None,
        "ldist_samp_lvl": None,
        "bc_seq": None,
        "filename": None,
    }

    for col, default in defaults.items():
        if col not in df.columns:
            df[col] = default

    # Accept legacy `run_num` on input, but normalize everything to `run_number`.
    if "run_num" in df.columns:
        df["run_number"] = df["run_number"].where(df["run_number"].notna(), df["run_num"])
        df = df.drop(columns=["run_num"])

    missing_run_id = df["run_id"].isna() | df["run_id"].astype(str).eq("")
    df.loc[missing_run_id, "run_id"] = df.loc[missing_run_id].apply(
        lambda row: make_run_id(
            row.get("run_number"),
            row.get("filename"),
            source_path,
        ),
        axis=1,
    )

    df["samp_date"] = df["samp_date"].apply(normalize_date)
    df["run_date"] = df["run_date"].apply(normalize_date)
    df["samp_group"] = df["samp_group"].map(normalize_label).astype("string")

    # Re-establish the numeric contract after CSV/XLSX round-trips.
    df["bc_count"] = pd.to_numeric(df["bc_count"], errors="coerce").fillna(0)
    df["proportion"] = pd.to_numeric(df["proportion"], errors="coerce").fillna(0.0)
    df["input"] = pd.to_numeric(df["input"], errors="coerce")
    df["x_idx"] = pd.to_numeric(df["x_idx"], errors="coerce")

    for col in ["above_cutoff", "multi_idx", "short_bc"]:
        df[col] = df[col].map(coerce_boolish)

    # Compile-mode identity is the composite biological sample id.
    missing_sample_id = df["sample_id"].isna() | df["sample_id"].astype(str).eq("")
    df.loc[missing_sample_id, "sample_id"] = df.loc[missing_sample_id].apply(
        lambda row: make_sample_id(
            row.get("samp_group"),
            row.get("samp_name"),
            row.get("samp_date"),
            row.get("run_number"),
            source_path,
        ),
        axis=1,
    )

    return df


def base_parent_name(s: pd.Series) -> pd.Series:
    """
    Strip mutation suffixes from parent annotations back to bare bc_name values.
    """
    s2 = s.astype("string")
    return s2.str.replace(r"_.*$", "", regex=True)

#%% Versions
"""
v3.0
- dropped matrix backend, refactored bk process
- added core_bc search
- changed param-passing to settings dict for flag_putative_parents() 
v2.2
- Removed higher-parent req for uniques
2026-05-26 - Refactored into bcparse package structure
"""
