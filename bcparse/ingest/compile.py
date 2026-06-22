# compile.py
"""
Name:       compile.py
Author:     CAG
Version:    2.0
Date:       2026/05/21
Refactored: 2026/05/26

Compile-mode ingest:
- parsed workbook and base CSV input normalization
- source deduplication
- core compiled-data QC
- `RunSeries` construction
"""

#%% Imports

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

import bcparse.qc as qc
from bcparse.containers.runinfo import RunInfo
from bcparse.containers.runseries import RunSeries

if TYPE_CHECKING:
    from bcparse.ingest.xlsx import ParsedAnalysis

#%% Compile functions

def read_base_csv(base_csv_path: str) -> pd.DataFrame:
    return pd.read_csv(base_csv_path).drop(
        columns=["ldist_group_lvl", "ldist_all"], errors="ignore"
    )

def build_runseries(
    *,
    parsed_files: list[ParsedAnalysis],
    settings: dict,
    base_df: pd.DataFrame | None = None,
) -> "RunSeries":
    """
    Build a deduplicated RunSeries from parsed analyses and optional existing compiled data.
    """
    compile_frames: list[pd.DataFrame] = []
    runinfo_by_run_id: dict[str, RunInfo] = {}

    if base_df is not None:
        csv_df = qc.normalize_long_df(
            base_df.drop(columns=["sample_id"], errors="ignore"),
            source_path=settings.get("base_csv_path"),
        )
        if not csv_df.empty:
            csv_df["__compile_source"] = "csv"
            csv_df["__compile_source_rank"] = 0
            csv_df["__compile_source_order"] = -1
            csv_df["__compile_source_mtime"] = float("-inf")
            compile_frames.append(csv_df)

    for run_order, parsed in enumerate(parsed_files):
        df = parsed.data_df.copy()
        if df.empty:
            continue

        run_ids = df["run_id"].dropna().astype(str).unique().tolist()
        if len(run_ids) != 1:
            raise ValueError(
                f"{parsed.filepath.name}: expected one run_id, found {run_ids}"
            )

        df["__compile_source"] = "analysis"
        df["__compile_source_rank"] = 1
        df["__compile_source_order"] = run_order
        df["__compile_source_mtime"] = parsed.filepath.stat().st_mtime

        compile_frames.append(df)
        runinfo_by_run_id[run_ids[0]] = parsed.runinfo

    if not compile_frames:
        combined_df = pd.DataFrame()
    else:
        print("- deduplicating parsed analyses and existing data...")
        combined_df = _dedup(pd.concat(compile_frames, ignore_index=True))

        print("- assembling final combined dataframe...")
        sort_cols = ["samp_group", "samp_name", "samp_date"]
        combined_df = (
            combined_df[sort_cols + [c for c in combined_df.columns if c not in sort_cols]]
            .sort_values(sort_cols, ascending=True, kind="stable")
            .reset_index(drop=True)
        )

        print("- running core QC...")
        combined_df = _core_qc(combined_df)

    n_source_runs = (
        combined_df["run_id"].dropna().astype(str).nunique()
        if not combined_df.empty else 0
    )

    run_series = RunSeries.from_long_df(
        combined_df,
        source="compiled",
        runinfo_by_run_id=runinfo_by_run_id,
    )

    base_csv_path = settings.get("base_csv_path")
    for run_id, seq_run in run_series.runs.items():
        seq_run.ensure_runinfo(filepath=base_csv_path)
        run_series.runinfo_by_run_id[run_id] = seq_run.runinfo

    run_series.series_meta.update({
        "out_prefix": settings.get("out_prefix"),
        "n_source_runs": n_source_runs,
    })
    run_series.data_df = combined_df.copy()

    return run_series

def _dedup(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate compile inputs by sample_id, preferring analysis over csv,
    then newer mtime, then later run order.
    """
    rank = ["__compile_source_rank", "__compile_source_mtime", "__compile_source_order"]

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
            source = df.loc[keep & df["sample_id"].eq(sid), "__compile_source"].iloc[0]
            print(f"  keeping: {sid} | source={source}")
            for source in group["__compile_source"].drop_duplicates():
                print(f"  dropped: {sid} | source={source}")

    deduped = df.merge(win, on=["sample_id", *rank], how="inner")

    meta_cols = [
        "sample_id", "idx_name", "run_id", "run_number", "run_name",
        "run_date", "samp_group", "samp_name", "samp_date", "input", "filename",
    ]
    blocks = deduped[[c for c in meta_cols if c in deduped.columns]].drop_duplicates()
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

def _core_qc(df: pd.DataFrame) -> pd.DataFrame:
    """Core QC steps applied to the combined compile dataframe."""
    df = df.copy()

    # Re-establish numeric columns after CSV/XLSX round-trips.
    df[["bc_count", "proportion", "input"]] = (
        df[["bc_count", "proportion", "input"]].apply(pd.to_numeric, errors="coerce").fillna(0)
    )

    df = qc.rank_uniques(df=df)
    df["bc_name"] = df["bc_name"].str.replace("SIVmac293M2", "SIVmac239M2", regex=False)
    df = qc.flag_shorts(df=df)

    df["ldist_samp_lvl"] = (
        df["ldist_samp_lvl"]
        .str.replace(r"^(1-)|(-1)$", "", regex=True)
        .str.replace(r"^(Unique)", "", regex=True)
        .replace("NA-NA", None)
    )

    df["multi_idx"] = (
        df.groupby(["run_number", "samp_group", "bc_name"])["samp_name"]
        .transform(lambda s: s.nunique() > 1)
    )

    # Annotate x_idx across the whole series to catch/refresh missing or wrong in old analyses.
    df = qc.annotate_x_idx(df)

    return df

#%% Versions
"""
v1.0 
- Initial version, accomodating new model/container structure
- Supercedes deprecated xlsx_compiler module and associated code in bcParse.py

v2.0
- Plain class replacing dataclass
- copy() and collapse_samples_in_place() removed; emit_qc works directly from sample_list
- Compile path collapsed from ~15 staticmethods to build_from_parsed_files + _dedup + _core_qc
- sample_id is the compile deduplication contract
2026-05-26 - Refactored into bcparse package structure
2026-05-26 - Ensured compiled runs receive RunInfo and updated compile readout
"""
