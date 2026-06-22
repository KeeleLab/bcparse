# flags.py
"""
Name:       flags.py
Author:     CAG
Version:    3.0
Date:       2026/03/20
Refactored: 2026/05/26

Consolidating qc functions for processing barcode data.
These are used across the parse- and compile-side build pipelines.
"""

#%% Imports
from __future__ import annotations

import numpy as np
import pandas as pd

#%% General

def get_prop(df: pd.DataFrame):
    """
    Note: written to run with single-index df via per_index_ops()
    """
    # Use internal copy
    df = df.copy()

    # calculate proportion
    df["proportion"] = df["bc_count"] / df["bc_count"].sum()

    return df


def flag_above_cutoff(df: pd.DataFrame):
    """
    Note: written to run with single-index df via per_index_ops()

    Assign above_cutoff as True for rows where:
        - bc_name does not contain 'Unique_' and proportion is greater than 1/input
        - bc_name contains 'Unique' and bc_count is greater than minimum named proportion * 100
    """
    # Use internal copy
    df = df.copy()

    # Initialize above_cutoff col
    df["above_cutoff"] = False

    # Pull given input
    input_val = df["input"].iloc[0]

    # Collect notes here
    notes = []

    # Handle infinite input gracefully
    if np.isinf(input_val) or np.isnan(input_val):
        named_cutoff = 1/10
        notes.append(f"Input value is {input_val}, used default cutoff {named_cutoff}")
    else:
        named_cutoff = 1 / input_val

    # Flag above_cutoff for non-unique barcodes
    df.loc[
        ~df["bc_name"].str.contains("Unique.")
        & (df["proportion"] > named_cutoff),
        "above_cutoff",
        ] = True

    # Get minimum named above_cutoff proportion
    min_named_p = df.loc[df["above_cutoff"], "proportion"].min()

    # Substitute nominal value if min_named_p is missing/nan
    if np.isnan(min_named_p):
        unique_cutoff = named_cutoff
        notes.append(f"No named above_cutoff: used prop > {unique_cutoff:.3g}")
    # Otherwise use 100X minimum named proportion
    else:
        unique_cutoff = 100 * min_named_p

    # Flag above_cutoff for unique barcodes
    df.loc[
        df["bc_name"].str.contains("Unique.")
        & (df["proportion"] > unique_cutoff),
        "above_cutoff",
        ] = True

    # Join all notes into a single string separated by semicolons
    note = "; ".join(notes) if notes else None

    return df, note


def rank_uniques(
    df: pd.DataFrame,
    bc_name_col: str = "bc_name",
    bc_seq_col: str = "bc_seq",
    bc_count_col: str = "bc_count"
    ) -> pd.DataFrame:
    """
    Rename 'Unique' barcodes by descending total count order,
    assigning names like 'Unique.1', 'Unique.2', etc.
    """
    for col in [bc_name_col, bc_seq_col, bc_count_col]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame.")

    # Use internal copy
    df = df.copy()

    # Build mapping dataframe with sum of counts per Unique barcode seq
    unique_mask = df[bc_name_col].astype(str).str.contains("Unique", na=False)
    unique_df = df.loc[unique_mask, [bc_seq_col, bc_count_col]].copy()
    mapping_df = (
        unique_df
        .groupby(bc_seq_col, as_index=False)
        .agg({bc_count_col: "sum"})
        .sort_values(by=bc_count_col, ascending=False)
        .reset_index(drop=True)
        )

    mapping_df[bc_name_col] = [f"Unique.{i+1}" for i in range(len(mapping_df))]

    # Create mapping dictionary from sequence to new unique name
    mapping = mapping_df.set_index(bc_seq_col)[bc_name_col].to_dict()

    # Apply mapping to 'Unique' barcode rows
    df.loc[unique_mask, bc_name_col] = df.loc[unique_mask, bc_seq_col].map(mapping)

    return df


def flag_shorts(
    df: pd.DataFrame,
    bc_seq_col: str = "bc_seq",
    motif: str = "ACTAGCATAA"
    ) -> pd.DataFrame:
    """
    Flags barcodes containing a "short" motif.
    Returns:
        DataFrame with a new boolean column 'short_bc'.
    """
    if bc_seq_col not in df.columns:
        raise ValueError(f"Column '{bc_seq_col}' not found in DataFrame.")

    # Use internal copy
    df = df.copy()

    df["short_bc"] = df[bc_seq_col].str.contains(motif, na=False)

    return df


def annotate_x_idx(df: pd.DataFrame) -> pd.DataFrame:
    """
    For shared barcodes, store each sample's share of that barcode's run total.
    Barcodes confined to one sample are left blank (NaN).
    """
    if df.empty:
        return df
    df = df.copy()
    group_cols = ["run_number", "bc_name"]
    n_samples = df.groupby(group_cols)["idx_name"].transform("nunique")
    bc_total  = df.groupby(group_cols)["bc_count"].transform("sum")
    df["x_idx"] = np.where(n_samples > 1, df["bc_count"] / bc_total.replace(0, np.nan), np.nan)
    return df

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
