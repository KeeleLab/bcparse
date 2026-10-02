# workbooks.py
"""
Name:       workbooks.py
Author:     CAG
Version:    2.1.0
Date:       20260623
"""

# %% Imports

from __future__ import annotations

import os
from datetime import datetime

import pandas as pd

import bcparse.qc as qc
from bcparse.containers.runseries import RunSeries
from bcparse.containers.seqrun import SeqRun

# %% Parse workbook


def write_parse_workbook(
    parse_run,
    *,
    settings: dict,
    runinfo_df: pd.DataFrame | None = None,
) -> None:
    """
    Write parse-mode csv/xlsx outputs from this parsed FASTQ run.
    """
    from bcparse.emit.views import samp_group_matrices, sidelong_tables

    out_path = settings["out_path"]
    legacy_format = settings["legacy_format"]
    filt_mat_ac = settings["filt_mat_ac"]
    collapse_to_parent = settings["collapse_to_parent"]
    seq_run = parse_run.seq_run
    runinfo_obj = parse_run.runinfo

    if seq_run is None:
        raise ValueError("Parse output requested before SeqRun was built.")

    # Initialize runinfo output table.
    if runinfo_df is None:
        runinfo_df = (
            pd.DataFrame() if runinfo_obj is None else runinfo_obj.copy_raw_df()
        )
    runinfo_df = runinfo_df.copy()
    runinfo_df["Date"] = runinfo_df["Date"].apply(qc.normalize_date)

    os.makedirs(out_path, exist_ok=True)

    # Write full long df for downstream re-use / inspection.
    concat_df = seq_run.df.copy()
    concat_df[["run_date", "samp_date"]] = concat_df[["run_date", "samp_date"]].apply(
        lambda col: col.apply(qc.normalize_date)
    )
    concat_df.to_csv(
        os.path.join(out_path, f"{seq_run.run_name}_concat.csv"), index=False
    )

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

    if not set(group_mats_dict.keys()) == set(groups) or not set(
        group_tabs_dict.keys()
    ) == set(groups):
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
    runinfo_df = _modify_runinfo(parse_run, runinfo_df, seq_run)

    workbook_kind = "Discovery" if settings.get("discover", False) else "Analysis"
    workbook_path = os.path.join(
        out_path,
        f"{seq_run.run_name}_{workbook_kind}.xlsx",
    )

    with pd.ExcelWriter(workbook_path) as writer:
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
        settings_df.to_excel(
            writer, sheet_name="Analysis settings", index=False, header=False
        )


def _modify_runinfo(
    parse_run,
    runinfo_df: pd.DataFrame,
    seq_run: SeqRun,
) -> pd.DataFrame:
    """
    Add parse summary columns and notes to the output runinfo table.
    """
    modified_runinfo = runinfo_df.copy()
    runinfo_obj = seq_run.runinfo
    index_key = (
        modified_runinfo["full_idx"]
        if runinfo_obj is None
        else runinfo_obj.index_key_for_output(modified_runinfo)
    )
    idx_notes = parse_run.notes_by_idx()

    modified_runinfo["n_reads"] = index_key.map(parse_run.total_reads_by_idx())
    modified_runinfo["n_bc"] = index_key.map(parse_run.total_barcodes_by_idx())
    modified_runinfo["n_bc_named"] = index_key.map(parse_run.named_barcodes_by_idx())
    modified_runinfo["n_reads_ac"] = index_key.map(
        parse_run.above_cutoff_reads_by_idx()
    )
    modified_runinfo["count-to-input_ratio"] = (
        modified_runinfo["n_reads"] / modified_runinfo["Input TOTAL PER BARCODE"]
    )

    if idx_notes:
        modified_runinfo["notes"] = index_key.map(idx_notes)

    return modified_runinfo


# %% Compile workbook


def write_compile_workbook(run_series: RunSeries, settings: dict) -> None:
    from bcparse.emit.views import (
        concat_series_runinfo,
        contam_report,
        group_mat_dict,
        sidelong_tables,
    )

    out_path, out_prefix = settings["out_path"], settings["out_prefix"]
    csv_path = os.path.join(out_path, f"{out_prefix}_all.csv")
    xlsx_path = os.path.join(out_path, f"{out_prefix}_above_cutoff.xlsx")
    os.makedirs(out_path, exist_ok=True)

    compile_df, df_ac, contam_check = emit_qc(run_series, settings=settings)
    compile_df_out = _drop_compile_private_cols(compile_df)
    df_ac_out = _drop_compile_private_cols(df_ac)

    date_cols = ["run_date", "samp_date"]
    compile_df_out[date_cols] = compile_df_out[date_cols].apply(
        lambda col: col.apply(qc.normalize_date)
    )
    df_ac_out[date_cols] = df_ac_out[date_cols].apply(
        lambda col: col.apply(qc.normalize_date)
    )

    compile_df_out.to_csv(csv_path, index=False)

    runinfo_df = concat_series_runinfo(run_series)

    sidelong_dict, _ = sidelong_tables(df_ac_out)
    _, mat_dict = group_mat_dict(df_ac_out)

    with pd.ExcelWriter(xlsx_path) as writer:
        if not runinfo_df.empty:
            runinfo_df.to_excel(writer, sheet_name="runinfo", index=False)

        df_ac_out.to_excel(writer, sheet_name="all_samples", index=False)

        if contam_check:
            full_mat, contam_table = contam_report(df_ac_out)
            full_mat.to_excel(writer, sheet_name="full_matrix", index=True)
            contam_table.to_excel(writer, sheet_name="contam_report", index=True)

        for samp_group in mat_dict:
            print(f"Writing {samp_group}...")
            sidelong_dict[samp_group].to_excel(
                writer,
                sheet_name=f"{samp_group} samples",
                index=False,
                header=False,
            )
            mat_dict[samp_group].to_excel(
                writer,
                sheet_name=f"{samp_group} matrix",
                index=True,
            )


def emit_qc(
    run_series: RunSeries, *, settings: dict
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    """Build QC-annotated compile outputs without mutating self."""
    if settings.get("collapse_to_parent", True):
        print("- collapsing children to parent within each sample...")
        compile_df = pd.concat(
            [
                sample.collapse_to_parent().df
                for run in run_series.run_list
                for sample in run.sample_list
            ],
            ignore_index=True,
        )
    else:
        compile_df = run_series.df.copy()

    if compile_df.empty:
        return pd.DataFrame(), pd.DataFrame(), False

    if any(
        settings.get(k, False)
        for k in ["ldist_samp_lvl", "ldist_group_lvl", "contam_check"]
    ):
        print("- running optional QC...")

    compile_df, contam_check = _apply_optional_qc(compile_df, settings)

    print("- building above-cutoff dataframe...")
    return compile_df, _build_above_cutoff_df(compile_df), contam_check


def _drop_compile_private_cols(df: pd.DataFrame) -> pd.DataFrame:
    private_cols = [col for col in df.columns if str(col).startswith("__compile_")]
    return df.drop(columns=private_cols, errors="ignore").copy()


# ====================
# QC HELPERS
# ====================


def _apply_optional_qc(df: pd.DataFrame, settings: dict) -> tuple[pd.DataFrame, bool]:
    if any(
        settings.get(k, False)
        for k in ["ldist_samp_lvl", "ldist_group_lvl", "contam_check"]
    ):
        print("\nRunning optional distance checks...")

    if settings.get("ldist_samp_lvl", False):
        df = _apply_ldist_qc(
            df,
            settings,
            groupby_cols=["samp_group", "samp_name", "samp_date"],
            output_col="ldist_samp_lvl",
            drop_existing=True,
        )

    if settings.get("ldist_group_lvl", False):
        df = _apply_ldist_qc(
            df,
            settings,
            groupby_cols=["samp_group"],
            aggr=True,
            output_col="ldist_group_lvl",
        )

    contam_check = bool(settings.get("contam_check", False))
    if contam_check:
        df = _apply_ldist_qc(
            df,
            settings,
            groupby_cols=None,
            aggr=True,
            output_col="ldist_all",
        )

    return df, contam_check


def _apply_ldist_qc(
    df: pd.DataFrame,
    settings: dict,
    groupby_cols: list[str] | None,
    aggr: bool = False,
    output_col: str = "putative_parent",
    drop_existing: bool = False,
) -> pd.DataFrame:
    if drop_existing and output_col in df.columns:
        df = df.drop(columns=[output_col])
    if groupby_cols:
        parts: list[pd.DataFrame] = []
        for group_key, group_df in df.groupby(groupby_cols, sort=True):
            set_name = " / ".join(map(str, group_key))
            parts.append(
                qc.flag_putative_parents(
                    set_name=set_name,
                    df=group_df,
                    settings=settings,
                    aggr=aggr,
                )
            )

        return pd.concat(parts, ignore_index=True).rename(
            columns={"putative_parent": output_col}
        )

    return qc.flag_putative_parents(
        set_name="Full Series", df=df, settings=settings, aggr=aggr
    ).rename(columns={"putative_parent": output_col})


def _build_above_cutoff_df(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "above_cutoff" not in df.columns:
        return df.iloc[0:0].copy()
    df_ac = df.loc[df["above_cutoff"]].copy()
    if df_ac.empty:
        return df_ac

    df_ac["bc_count"] = pd.to_numeric(df_ac["bc_count"], errors="coerce").fillna(0)
    df_ac["__denom"] = df_ac.groupby("sample_id", dropna=False)["bc_count"].transform(
        "sum"
    )
    df_ac["proportion"] = (df_ac["bc_count"] / df_ac["__denom"]).fillna(0)

    return df_ac.drop(columns=["__denom"])


# %% Versions
"""
v2.1.0 20260623
 - accomodating source/sample_id meta for views.py update

v2.0.0 20260526
 - Refactored into bcparse package structure
 - Updated compile optional-QC console reporting
"""
