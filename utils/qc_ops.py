# qc_ops.py
"""
Names:      qc_ops.py
Author:     CAG
Version:    3.0
Date:       2026/03/20

Consolidating qc fucntions for processing barcode data.
These are used in parse_countdata.py and qc_ops.py
"""

#%% Imports
from __future__ import annotations
from typing import Sequence, Literal, Dict, List, Optional, Callable, Tuple
import numpy as np
import pandas as pd
from dataclasses import dataclass
from rapidfuzz.distance import Levenshtein
from utils.bk_tree import build_ref_index, query

#%% General

def normalize_date(date_val) -> Optional[str]:
    """
    Coerce a value to a normalized date string (YYYY-MM-DD) if possible.
    """
    try:
        dt = pd.to_datetime(date_val, errors="coerce")
        if pd.isna(dt):
            return None
        return dt.strftime("%Y-%m-%d")
    except Exception:
        return None


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

    # Build mapping dataframe with sum of counts per barcode seq
    mapping_df = (
        df[df[bc_name_col].str.contains("Unique", na=False)]
        .groupby(bc_seq_col, as_index=False)[bc_count_col]
        .sum()
        .sort_values(bc_count_col, ascending=False)
        .reset_index(drop=True)
        )

    mapping_df[bc_name_col] = [f"Unique.{i+1}" for i in range(len(mapping_df))]

    # Create mapping dictionary from sequence to new unique name
    mapping = mapping_df.set_index(bc_seq_col)[bc_name_col].to_dict()

    # Apply mapping to 'Unique' barcode rows
    mask = df[bc_name_col].str.contains("Unique", na=False)
    df.loc[mask, bc_name_col] = df.loc[mask, bc_seq_col].map(mapping)

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


#%% Parent ID 

def flag_putative_parents(
    settings: dict,
    set_name: str,
    df: 'pd.DataFrame',
    bc_name_col: str = "bc_name",
    bc_seq_col: str = "bc_seq",
    bc_count_col: str = "bc_count",
    *,
    # optional args for multi-sample check
    aggr: bool = False,
    proportion_col: str = "proportion",
    input_col: str = "input",
    ) -> 'pd.DataFrame':
    """
    Annotate df with 'putative_parent' column indicating lineage origin.
    """

    # Pull vars from settings
    dist = settings.get("dist_threshold")
    core_bc = settings.get("core_bc")

    # Define working df
    df_c = df.copy()

    # Optional preprocessing for `aggr` runs
    if aggr:

        # Replace inf or NaN with 10
        df_c[input_col] = df_c[input_col].replace([np.inf, np.nan], 10)

        # Calculate adjusted counts, round up and set as int
        df_c["adj_count"] = np.ceil(df_c[proportion_col] * df_c[input_col]).astype(int)

        # Aggregate adjusted counts by barcode
        df_c = (
            df_c.groupby([bc_name_col, bc_seq_col], as_index=False)
            .agg({'adj_count': 'sum'})
            .rename(columns={'adj_count': bc_count_col})
            )

    # Subset named and unique
    df_n = df_c[~df_c[bc_name_col].str.startswith("Unique", na=False)]
    df_u = df_c[df_c[bc_name_col].str.startswith("Unique", na=False)]

    # Get numpy arrays for vectorized / indexed access
    # Named
    n_seqs   = df_n[bc_seq_col].astype(str).to_numpy()
    n_counts = df_n[bc_count_col].astype(int).to_numpy()
    n_names  = df_n[bc_name_col].astype(object).to_numpy()
    # Unique
    u_seqs   = df_u[bc_seq_col].astype(str).to_numpy()
    u_counts = df_u[bc_count_col].astype(int).to_numpy()
    u_names  = df_u[bc_name_col].astype(object).to_numpy()

    # Build ref tree on named
    ref_tree = build_ref_index(
        ref_ids = np.arange(len(n_names), dtype=int).tolist(),
        ref_seqs = n_seqs.tolist(), 
        distance = "Levenshtein", 
        verify = False
    )

    # Query named vs. named
    res_named = query(
        tree = ref_tree,
        query_ids = list(n_names),
        query_seqs = list(n_seqs),
        mode="single",
        scope="radius",
        max_radius=int(dist),
        returns="pairs",  # [df_q_idx[(df_r_idx, dist), ...]]
        ) if not df_n.empty else []

    # Pass to parent select (req p>c, suffix shows subs)
    named_parents = _select_parents(
        q_names=n_names,
        q_seqs=n_seqs,
        q_counts=n_counts,
        q_hits=res_named,
        r_names=n_names,
        r_seqs=n_seqs,
        r_counts=n_counts,
        req_parent_higher=True,
        suffix_fn=_parent_suffix_sub,
        )
 
    # Query unique vs. named
    res_unique = query(
        tree = ref_tree,
        query_ids = list(u_names),
        query_seqs = list(u_seqs),
        mode="single",
        scope="radius",
        max_radius=int(dist),
        returns="pairs",  # [df_q_idx[(df_r_idx, dist), ...]]
        ) if not df_u.empty else []

    # Pass to parent select (!req p>c, suffix shows subs)
    unique_parents = _select_parents(
        q_names=u_names,
        q_seqs=u_seqs,
        q_counts=u_counts,
        q_hits=res_unique,
        r_names=n_names,
        r_seqs=n_seqs,
        r_counts=n_counts,
        req_parent_higher=False,
        suffix_fn=_parent_suffix_sub,
        )

    # Build parent mapping dict
    all_parents = {}
    all_parents.update(named_parents)
    all_parents.update(unique_parents)

    # Optional: check remaining uniques via sliding-window core bc search
    if core_bc is not None: 

        # Pull idx vals for unresolved uniques and search
        unresolved_idx = np.flatnonzero(np.array([len(hits) == 0 for hits in res_unique], dtype=bool))
        df_u2 = df_u.iloc[unresolved_idx].reset_index(drop=True)

        if len(unresolved_idx) > 0:
            res_unique_core = _core_bc_search(
                unresolved_df=df_u2,
                ref_seqs=n_seqs,  # search against named barcodes
                core_bc=core_bc,
                bc_seq_col=bc_seq_col,
                )

            u2_seqs   = df_u2[bc_seq_col].astype(str).to_numpy()
            u2_counts = df_u2[bc_count_col].astype(int).to_numpy()
            u2_names  = df_u2[bc_name_col].astype(object).to_numpy()

            # Pass to parent selection (req p>c, suffix shows match locus)
            unique_parents_core = _select_parents(
                q_names=u2_names,
                q_seqs=u2_seqs,
                q_counts=u2_counts,
                q_hits=res_unique_core,
                r_names=n_names,
                r_seqs=n_seqs,
                r_counts=n_counts,
                req_parent_higher=False,
                suffix_fn=_parent_suffix_full,
                )
            
            if unique_parents_core:
                all_parents.update(unique_parents_core)

    df["putative_parent"] = df[bc_name_col].map(all_parents)

    return df


def _core_bc_search(
    *,
    unresolved_df: pd.DataFrame,
    ref_seqs: np.ndarray,
    core_bc: tuple[int, int],
    bc_seq_col: str = "bc_seq",
    ) -> list[list[tuple[int, int]]]: # [df_q_idx[(df_r_idx, dist), ...]]
    """
    Sliding-window fallback search over the core barcode region.
    """
    core_start, core_len = core_bc
    core_end = core_start + core_len

    # Build reference index on the core_bc slice
    ref_seqs_core = np.array([s[core_start:core_end] for s in ref_seqs], dtype=object)

    core_tree = build_ref_index(
        ref_ids=np.arange(len(ref_seqs_core), dtype=int).tolist(),
        ref_seqs=ref_seqs_core.tolist(),
        distance="Levenshtein",
        verify=False,
    )

    # Init res list
    q_res: list[list[tuple[int, int]]] = []

    # sliding-window search
    for seq in unresolved_df[bc_seq_col].astype(str).to_numpy():
        if len(seq) < core_len:
            q_res.append([])
            continue

        # Per query, define window slice lists
        win_q_ids = [f"0:{start}" for start in range(0, len(seq) - core_len + 1)]
        win_q_seqs = [seq[start : start + core_len] for start in range(0, len(seq) - core_len + 1)]

        # Query that list
        win_hits = query(
            tree=core_tree,
            query_ids=win_q_ids,
            query_seqs=win_q_seqs,
            mode="single",
            scope="radius",
            max_radius=0,  # Require exact match for core searches
            returns="pairs",
        )

        # Collape to single list per query
        seen: set[int] = set()
        core_hits: list[tuple[int, int]] = []

        for pairs in win_hits:
            for ref_idx, dist in pairs:
                if ref_idx not in seen:
                    seen.add(ref_idx)
                    core_hits.append((ref_idx, dist))

        q_res.append(core_hits)

    return q_res


def _select_parents(
    *,
    q_names: np.ndarray,
    q_seqs: np.ndarray,
    q_counts: np.ndarray,
    q_hits: list[list[tuple[int, int]]],
    r_names: np.ndarray,
    r_seqs: np.ndarray,
    r_counts: np.ndarray,
    req_parent_higher: bool,
    suffix_fn: Optional[Callable[[str, str], str]] = None,
    ) -> dict[str, str]:
    """
    Resolve to one parent per query.
    """
    if suffix_fn is None:
        suffix_fn = _parent_suffix_full

    r_order = np.arange(len(r_names), dtype=int)
    results: dict[str, str] = {}

    for q_name, q_seq, q_count, hits in zip(q_names, q_seqs, q_counts, q_hits):
        
        # skip non-hits
        if not hits:
            continue

        # get reference index
        cand_ref_idx = np.fromiter((rid for rid, _ in hits), dtype=int)

        # drop self-matches/continue
        not_self = r_names[cand_ref_idx] != q_name
        if not np.any(not_self):
            continue

        cand_ref_idx = cand_ref_idx[not_self]

        if req_parent_higher:
            # Select highest parent w/ count > query
            # Filter out parents < query...
            count_diffs = r_counts[cand_ref_idx] - int(q_count)
            keep = count_diffs > 0
            if not np.any(keep):
                continue
            # Choose highest count, tiebreak by list order
            cand_ref_idx = cand_ref_idx[keep]
            count_diffs = count_diffs[keep]
            ref_rank = r_order[cand_ref_idx]
            pick = np.lexsort((ref_rank, -count_diffs))[0]

        else:
            # Just select highest available parent
            ref_counts = r_counts[cand_ref_idx]
            ref_rank = r_order[cand_ref_idx]
            pick = np.lexsort((ref_rank, -ref_counts))[0]

        best_ref_idx = int(cand_ref_idx[pick])
        parent_name = str(r_names[best_ref_idx])
        suffix = suffix_fn(r_seqs[best_ref_idx], q_seq)

        results[str(q_name)] = f"{parent_name}_{suffix}" if suffix else parent_name

    return results


def _parent_suffix_full(parent_seq: str, child_seq: str) -> str:
    """
    Return a full CIGAR-like suffix for parent -> child.
    """
    ops = Levenshtein.editops(parent_seq, child_seq)
    if not ops:
        return ""
    
    # init local vars
    parts: list[str] = []
    cur_op = None
    cur_n = 0
    
    # Helper: Run-length encode operations, flush when the op changes.
    def emit(op: str) -> None:
        nonlocal cur_op, cur_n
        if op == cur_op:
            cur_n += 1
            return
        if cur_op is not None:
            parts.append(f"{cur_n}{cur_op}")
        cur_op, cur_n = op, 1

    # Map Levenshtein ops → (CIGAR op, parent step, child step)
    step = {
        "replace": ("X", 1, 1),
        "insert": ("I", 0, 1),
        "delete": ("D", 1, 0),
        }

    i = j = 0 # pointers into parent_seq / child_seq
    
    for op in ops:
        # Fill exact matches up to next edit boundary.
        # Invariant: (i, j) tracks next unprocessed positions.
        while i < op.src_pos and j < op.dest_pos:
            emit("=")
            i += 1
            j += 1
        # Apply the edit op and advance indices
        cigar_op, di, dj = step[op.tag]
        emit(cigar_op)
        i += di
        j += dj
    
    # Emit any trailing matches after the final edit.
    while i < len(parent_seq) and j < len(child_seq):
        emit("=")
        i += 1
        j += 1
    
    # Flush final run
    if cur_op is not None:
        parts.append(f"{cur_n}{cur_op}")
    return "".join(parts)


def _parent_suffix_sub(parent_seq: str, child_seq: str) -> str:
    """
    Return compact substitution suffix for parent->child, e.g. "5C_12A".
    """
    ops = Levenshtein.editops(parent_seq, child_seq)
    subs = []
    for op in ops:
        if op.tag == "replace":
            pos = op.src_pos + 1             # 1-based parent position
            alt = child_seq[op.dest_pos]     # child's base at that position
            subs.append(f"{pos}{alt}")
    return "-".join(subs) if subs else ""


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
    def base_parent_name(s: pd.Series) -> pd.Series:
        # coerce to pandas string dtype (preserves missing as <NA>)
        s2 = s.astype("string")
        # safe replace; .str is valid on string dtype
        return s2.str.replace(r"_.*$", "", regex=True)

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
    Collapse all children to their parent rows.
    """
    if df["putative_parent"].isna().all():
        return df, None

    df = df.copy()

    def base_parent_name(s: pd.Series) -> pd.Series:
        s2 = s.astype("string")
        return s2.str.replace(r"_.*$", "", regex=True)

    parent_base = base_parent_name(df["putative_parent"])
    has_parent = df["putative_parent"].notna()

    child_rows = df.loc[has_parent].copy()
    if child_rows.empty:
        return df, None

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

#%% Versiont
"""
v3.0
- dropped matrix backend, refactored bk process
- added core_bc search
- changed param-passing to settings dict for flag_putative_parents() 
v2.2
- Removed higher-parent req for uniques
"""