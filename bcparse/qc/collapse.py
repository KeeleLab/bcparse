# collapse.py
"""
Name:       collapse.py
Author:     CAG
Version:    3.0
Date:       2026/03/20
Refactored: 2026/05/26

Consolidating qc functions for processing barcode data.
These are used across the parse- and compile-side build pipelines.
"""

#%% Imports
from __future__ import annotations
import pandas as pd
from bcparse.qc.identity import base_parent_name
from bcparse.qc.flags import get_prop, flag_above_cutoff

#%% Collapsing functions

def collapse_ambig_children(df: pd.DataFrame):
    """
    Collapse children with a single 'N' to valid parent rows, and drop multi-'N'.
    
    Assumes:
    - 'bc_seq' (string with possible 'N's)
    - 'bc_name' (parent name without suffix)
    - 'putative_parent' (now bc_name + '-#' suffix)
    - 'above_cutoff' (bool)
    - 'bc_count' (numeric)
    Returns: (df, note)
    """
    # quick exit if no Ns exist at all
    if not df["bc_seq"].str.contains("N").any():
        return df, None

    note = "ambig_children collapsed"
    df = df.copy()

    # helper: strip a single trailing -<suffix> from the putative_parent
    # Drop rows where 'bc_seq' contains more than one 'N'
    df = df[~(df["bc_seq"].str.count("N") > 1)]

    # Identify ambiguous children
    parent_base = base_parent_name(df["putative_parent"].astype("string"))
    eligible_parents = df.loc[df["above_cutoff"], "bc_name"]

    ambig_mask = (
        df["putative_parent"].notna()         # has named putative parent
        & df["bc_seq"].str.contains("N")      # contains a single N (multi was dropped)
        & (~df["above_cutoff"])               # is below cutoff
        & parent_base.isin(eligible_parents)  # parent is above cutoff
        )

    ambig_children = df.loc[ambig_mask].copy()

    # Early exit if nothing to do
    if ambig_children.empty:
        return df, None

    # Remove these plus all remaining rows that contain any N
    df = df[~df["bc_seq"].str.contains("N")]

    # Sum ambiguous children counts by *base* parent name (bc_name without suffix)
    add_counts = ambig_children["bc_count"].groupby(base_parent_name(ambig_children["putative_parent"])).sum()

    # Add those counts to the corresponding parent rows
    df.loc[df["bc_name"].isin(add_counts.index), "bc_count"] += (
        df["bc_name"].map(add_counts).fillna(0)
        )

    # Recompute proportions
    df = get_prop(df)

    return df, note


def collapse_to_parent(df: pd.DataFrame):
    """
    Collapse all children to their parent rows within the scope of the
    dataframe provided.

    In other words, this function does not define grouping boundaries itself;
    it assumes the caller has already subset the dataframe to the intended
    collapse level (e.g. one SeqSamp, one idx_name, one group aggregate, etc).
    """
    if "putative_parent" in df.columns and df["putative_parent"].notna().any():
        parent_col = "putative_parent"
    elif "ldist_samp_lvl" in df.columns and df["ldist_samp_lvl"].notna().any():
        parent_col = "ldist_samp_lvl"
    else:
        return df

    df = df.copy()

    parent_base = base_parent_name(df[parent_col])
    has_parent = df[parent_col].notna()

    child_rows = df.loc[has_parent].copy()
    if child_rows.empty:
        return df

    # keep all non-child rows, including the parent rows themselves
    df = df.copy()

    add_counts = child_rows["bc_count"].groupby(parent_base.loc[has_parent]).sum()

    parent_mask = df["bc_name"].isin(add_counts.index)
    df.loc[parent_mask, "bc_count"] += df.loc[parent_mask, "bc_name"].map(add_counts).fillna(0)

    # drop the child rows after counts have been merged into parents
    df = df.loc[~has_parent].copy()

    # recycle these from above
    df = get_prop(df)
    df, _ = flag_above_cutoff(df)

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
