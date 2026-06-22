# parents.py
"""
Name:       parents.py
Author:     CAG
Version:    3.0
Date:       2026/03/20
Refactored: 2026/05/26

Consolidating qc functions for processing barcode data.
These are used across the parse- and compile-side build pipelines.
"""

#%% Imports
from __future__ import annotations
from typing import Callable, Optional
import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from bcparse.lib.bk_tree import build_ref_index, query

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
    print(f"{set_name}", flush=True)
    print(f" - building distance index", flush=True)
    ref_tree = build_ref_index(
        ref_ids = np.arange(len(n_names), dtype=int).tolist(),
        ref_seqs = n_seqs.tolist(), 
        distance = "Levenshtein", 
        verify = False
    )

    # Query named vs. named
    print(f" - querying named barcodes", flush=True)
    res_named = query(
        tree = ref_tree,
        query_ids = list(n_names),
        query_seqs = list(n_seqs),
        mode="single",
        scope="radius",
        max_radius=int(dist),
        returns="pairs",  # [df_q_idx[(df_r_idx, dist), ...]]
        ) if not df_n.empty else []

    # Pass to parent select (req p>c, suffix shows ldist ops)
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
    print(f" - querying unique barcodes", flush=True)
    res_unique = query(
        tree = ref_tree,
        query_ids = list(u_names),
        query_seqs = list(u_seqs),
        mode="single",
        scope="radius",
        max_radius=int(dist),
        returns="pairs",  # [df_q_idx[(df_r_idx, dist), ...]]
        ) if not df_u.empty else []

    # Pass to parent select (!req p>c, suffix shows ldist ops)
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

            # Pass to parent selection (req p>c, tag as core fallback hit)
            unique_parents_core = _select_parents(
                q_names=u2_names,
                q_seqs=u2_seqs,
                q_counts=u2_counts,
                q_hits=res_unique_core,
                r_names=n_names,
                r_seqs=n_seqs,
                r_counts=n_counts,
                req_parent_higher=False,
                suffix_fn=_parent_suffix_core,
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
            max_radius=1,  # Allow one edit in sliding-window core searches
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
        return f"{len(parent_seq)}="
    
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


def _parent_suffix_core(parent_seq: str, child_seq: str) -> str:
    """
    Return the fixed tag used for core-barcode fallback hits.
    """
    return "core"


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

#%% Versions
"""
v3.1
- Updatated console reporting
v3.0
- dropped matrix backend, refactored bk process
- added core_bc search
- changed param-passing to settings dict for flag_putative_parents() 
v2.2
- Removed higher-parent req for uniques
"""
