#!/usr/bin/env python3
# bk_tree.py

"""
Name:       bcparse
Author:     CAG
Version:    1.0.1
Date:       20240605

BK-tree: Approximate string matching under edit or Hamming distance
===================================================================
A lightweight, self-contained implementation of a BK-tree index for
rapid nearest-neighbor search over strings using Levenshtein or
Hamming distance (RapidFuzz backend). Designed for integration in
data pipelines or ad-hoc use from the Python CLI.


Primary entry points
--------------------
- build_ref_index(ref_ids, ref_seqs, distance="Levenshtein", verify=False)
    Build a BK-tree index once. Returns a BKTree object.

- query(tree, query_ids, query_seqs, *, mode="single", scope="min", ...)
    Perform flexible nearest-neighbor queries. Supports multiple
    result shapes, tie-breaking modes, and radius-based searches.


Common parameters
-----------------
scope:
    "min"      return only neighbors at the minimal distance
    "radius"   all neighbors with distance <= max_radius
    "exact"    all neighbors with distance == max_radius

mode:
    "one"      single query, single return
    "single"   multiple queries, one result each
    "all"      exhaustive pairwise (for edge/graph building)

returns:
    "parallel" aligned lists ([ids], [dists])
    "pairs"    list of (id, dist) tuples
    "map"      {id: dist} dictionary
    "edges"    list of (query_id, ref_id, dist) triples (scope ≠ "min")

tie:
    "first" | "min_id" | "all" — controls behavior for equal minima


mode     scope     returns     shape
-----    -------   --------    ---------------------------------------------
one      min       (ignored)   (best_id or [ids if tie='all'], dist_or_None)
one      radius    pairs       [(id, dist)]
one      radius    parallel    ([ids], [dists])
one      radius    map         {id: dist}
single   min       -           ([best_id?], [best_dist?])  # aligned with queries
single   radius    pairs       [[(id, dist)], ...]         # per-query
all      radius    edges       [(query_id, id, dist)]      # flat edge list


Typical CLI or REPL usage
-------------------------
from bk_tree import build_ref_index, query
refs = ["apple", "apply", "maple", "snapple"]
rids = ["r1", "r2", "r3", "r4"]
tree = build_ref_index(rids, refs)
query(tree, ["q1"], ["appl"], mode="one", scope="min")
('r1', 1)


Notes
-----
* Fully deterministic: multi-hit results are sorted by (distance, id).
* No command-line interface required — import or `python -c` one-liners.
* Designed for composability in larger QC, barcode, or clustering pipelines.
"""

import random
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

from rapidfuzz import distance as rfdist

# ----------
# Result types
# ----------

# Added to support type checking and clarity in return types for the query function.
# The alternative would be to have seperate functions for each return type.

MinOneResult = tuple[Any, int | None]
MinAllResult = tuple[list[Any], int | None]

OneParallelResult = tuple[list[Any], list[int]]
OnePairsResult = list[tuple[Any, int]]
OneMapResult = dict[Any, int]

BatchMinResult = tuple[list[Any | None], list[int | None]]
BatchParallelResult = tuple[list[list[Any]], list[list[int]]]
BatchPairsResult = list[list[tuple[Any, int]]]
BatchMapResult = list[dict[Any, int]]
EdgeResult = list[tuple[Any, Any, int]]

NearestResult = (
    MinAllResult
    | MinOneResult
    | OneParallelResult
    | OnePairsResult
    | OneMapResult
)

QueryResult = (
    MinOneResult
    | MinAllResult
    | OneParallelResult
    | OnePairsResult
    | OneMapResult
    | BatchMinResult
    | BatchParallelResult
    | BatchPairsResult
    | BatchMapResult
    | EdgeResult
)

# ----------
# Distance
# ----------

class DistanceSpec:
    """
    A self-contained distance "object":
      - Construct from a callable or a string (e.g., 'Levenshtein','lev','ham','Indel','OSA',...)
      - Callable via __call__(a,b) -> int
      - Provides probe() / assert_ok() for quick metric checks.
    """

    def __init__(self, distance: Union[str, Callable[[str, str], Any]]):
        self._name: Optional[str] = None
        self._fn: Callable[[str, str], int] = self._resolve(distance)

    # ---- public API ----
    def __call__(self, a: str, b: str) -> int:
        # Always return a plain int
        return int(self._fn(a, b))

    @property
    def name(self) -> str:
        return self._name or "custom"

    # ---- helpers / resolution ----
    @staticmethod
    def _normalize_name(s: str) -> str:
        return s.replace("_", "").replace("-", "").lower()

    @staticmethod
    def _wrap_hamming(fn: Callable[[str, str], Any]) -> Callable[[str, str], int]:
        def _h(a: str, b: str) -> int:
            if len(a) != len(b):
                raise ValueError("Hamming distance requires equal-length strings")
            return int(fn(a, b))

        return _h

    def _resolve(
        self, distance: Union[str, Callable[[str, str], Any]]
    ) -> Callable[[str, str], int]:
        # Direct callable
        if callable(distance):
            # Light sanity: must be int-convertible on trivial input
            d0 = distance("a", "a")
            _ = int(d0)  # may raise -> good
            self._name = getattr(distance, "__name__", "custom")
            return lambda a, b: int(distance(a, b))

        # String → rapidfuzz distance
        if isinstance(distance, str):
            name = self._normalize_name(distance)
            self._name = distance

            # Common aliases
            if name in ("lev", "levenshtein", "edit"):
                return lambda a, b: int(rfdist.Levenshtein.distance(a, b))
            if name in ("ham", "hamming"):
                return self._wrap_hamming(rfdist.Hamming.distance)

            # Case-insensitive lookup under rapidfuzz.distance
            candidates = {
                k: getattr(rfdist, k) for k in dir(rfdist) if not k.startswith("_")
            }
            match = next(
                (v for k, v in candidates.items() if self._normalize_name(k) == name),
                None,
            )
            if match is None or not hasattr(match, "distance"):
                raise ValueError(
                    f"Unknown distance '{distance}'. "
                    "Use a rapidfuzz distance name (e.g., 'Levenshtein','Hamming','Indel','OSA','LCSseq', ...) "
                    "or provide a callable."
                )
            # Special-case Hamming-like
            if self._normalize_name(match.__name__) == "hamming":
                return self._wrap_hamming(match.distance)

            # Generic: ensure int-convertible
            d0 = match.distance("a", "a")
            _ = int(d0)
            return lambda a, b: int(match.distance(a, b))

        raise ValueError(
            "distance must be a callable or a string name of a rapidfuzz distance"
        )

    # ---- metric checks ----
    def probe(
        self,
        seqs: Sequence[str],
        *,
        n_pairs: int = 200,
        n_triples: int = 300,
        rng_seed: int = 42,
    ) -> Dict[str, int]:
        """
        Sample-based metric check.
        Returns dict: {'nonint': x, 'identity': y, 'symmetry': z, 'triangle': w}
        """
        rng = random.Random(rng_seed)
        N = len(seqs)
        if N < 2:
            return {"nonint": 0, "identity": 0, "symmetry": 0, "triangle": 0}

        vio = {"nonint": 0, "identity": 0, "symmetry": 0, "triangle": 0}
        # Identity/type
        for _ in range(min(n_pairs, N)):
            a = seqs[rng.randrange(N)]
            try:
                d_aa = int(self(a, a))
            except Exception:
                vio["nonint"] += 1
                continue
            if d_aa != 0:
                vio["identity"] += 1

        # Symmetry
        for _ in range(min(n_pairs, N * (N - 1))):
            i, j = rng.randrange(N), rng.randrange(N)
            if i == j:
                continue
            a, b = seqs[i], seqs[j]
            try:
                dab = int(self(a, b))
                dba = int(self(b, a))
            except Exception:
                vio["nonint"] += 1
                continue
            if dab != dba:
                vio["symmetry"] += 1

        # Triangle
        for _ in range(n_triples):
            i, j, k = rng.randrange(N), rng.randrange(N), rng.randrange(N)
            a, b, c = seqs[i], seqs[j], seqs[k]
            try:
                dab = int(self(a, b))
                dbc = int(self(b, c))
                dac = int(self(a, c))
            except Exception:
                vio["nonint"] += 1
                continue
            if dac > dab + dbc:
                vio["triangle"] += 1

        return vio

    def assert_ok(
        self,
        seqs: Sequence[str],
        *,
        max_nonint: int = 0,
        max_identity: int = 0,
        max_symmetry: int = 0,
        max_triangle: int = 0,
        n_pairs: int = 200,
        n_triples: int = 300,
        rng_seed: int = 42,
    ) -> None:
        vio = self.probe(seqs, n_pairs=n_pairs, n_triples=n_triples, rng_seed=rng_seed)
        if (
            vio["nonint"] > max_nonint
            or vio["identity"] > max_identity
            or vio["symmetry"] > max_symmetry
            or vio["triangle"] > max_triangle
        ):
            raise ValueError(
                "Distance metric check failed: "
                f"{vio} (allowed: "
                f"nonint≤{max_nonint}, identity≤{max_identity}, symmetry≤{max_symmetry}, triangle≤{max_triangle})"
            )


# --------
# BK-tree
# --------


class BKNode:
    __slots__ = ("key", "id", "children")

    def __init__(self, key: str, id: Any):
        self.key: str = key
        self.id: Any = id
        # children keyed by integer edit distance from this node's key
        self.children: Dict[int, "BKNode"] = {}


class BKTree:
    def __init__(
        self,
        distance: Union[str, Callable[[str, str], Any], DistanceSpec] = "Levenshtein",
    ):
        self.root: Optional[BKNode] = None
        self.dist = (
            distance if isinstance(distance, DistanceSpec) else DistanceSpec(distance)
        )

    # --- Dunder methods for introspection and REPL friendliness ---
    def __len__(self) -> int:
        """
        Return the number of nodes (reference sequences) in the tree.
        Allows using len(tree) for quick inspection.
        """
        if self.root is None:
            return 0
        n = 0
        stack = [self.root]
        while stack:
            node = stack.pop()
            n += 1
            stack.extend(node.children.values())
        return n

    def __bool__(self) -> bool:
        """
        Return True if the tree has any nodes.
        Enables 'if tree:' semantics in pipelines.
        """
        return self.root is not None

    def __repr__(self) -> str:
        """
        Return a concise, human-readable summary of the tree.
        Shows size and distance metric when printed or inspected.
        Example:
            >>> tree
            <BKTree size=10234 distance='Levenshtein'>
        """
        name = getattr(self.dist, "name", "custom")
        return f"<BKTree size={len(self)} distance={name!r}>"

    # ---- construct ----
    def insert(self, key: str, id: Any):
        """Insert (sequence, id) into the BK-tree."""
        if self.root is None:
            self.root = BKNode(key, id)
            return
        node = self.root
        while True:
            d = self.dist(key, node.key)
            child = node.children.get(d)
            if child is None:
                node.children[d] = BKNode(key, id)
                return
            node = child

    def build(self, keys: Iterable[str], ids: Iterable[Any]):
        """Bulk insert; keys and ids must be aligned iterables of equal length."""
        for k, i in zip(keys, ids):
            self.insert(k, i)

    # ---- query ----
    def nearest(
        self,
        query_seq: str,
        query_id: Any = None,
        *,
        max_radius: Optional[int] = None,
        exclude_same_id: bool = True,
        tie: str = "all",  # "all" | "first" | "min_id"
        scope: str = "min",  # "min" | "radius" | "exact"
        returns: str = "parallel",  # "parallel" | "pairs" | "map" | "single"
    ) -> NearestResult:
        """
        Search the BK-tree for neighbors of query_seq.

        Parameters
        ----------
        query_seq : str
            Query sequence.
        query_id : Any, optional
            Label to exclude (prevents self-match when querying against same ref set).
        max_radius : int, optional
            Required for scope="radius" and scope="exact".
        exclude_same_id : bool
            If True, suppress candidate whose id == query_id.
        tie : {"all","first","min_id"}
            Used only when scope="min". Controls how to collapse equal-best ties.
        scope : {"min","radius","exact"}
            - "min":    return only neighbors at the minimal distance (<= max_radius if provided)
            - "radius": return ALL neighbors with distance <= max_radius   (requires max_radius)
            - "exact":  return ALL neighbors with distance == max_radius   (requires max_radius)
        returns : {"parallel","pairs","map","single"}
            For scope="radius"/"exact":
              - "parallel": returns two aligned lists (ids, dists)
              - "pairs":    returns list of (id, dist) tuples
              - "map":      returns dict {id: dist}  (dedup by id)
            For scope="min":
              - "single":   returns (one_id, dist) with tie policy
              - "parallel"/"pairs"/"map" fall back to the natural shape for scope="min":
                   tie="all" -> (list_ids, dist)
                   tie!="all"-> (one_id, dist)

        Returns
        -------
        See return type annotations above.
        """
        scope = scope.lower()
        tie = tie.lower()
        returns = returns.lower()

        if scope not in ("min", "radius", "exact"):
            raise ValueError('scope must be "min", "radius", or "exact"')

        if scope in ("radius", "exact") and max_radius is None:
            raise ValueError(f'scope="{scope}" requires max_radius')

        # Returning an empty query result when the tree has no root
        if self.root is None:
            if scope in ("radius", "exact") or tie == "all":
                # empty collection
                return (
                    ([], None)
                    if scope == "min"
                    else ([], [])
                    if returns == "parallel"
                    else []
                    if returns == "pairs"
                    else {}
                )
            return (None, None)

        # Collectors for scope="min"
        best_ids: List[Any] = []
        best_d: Optional[int] = None

        # Collectors for multi-hit scopes (radius/exact)
        ids_multi: List[Any] = []
        dists_multi: List[int] = []
        pairs_multi: List[Tuple[Any, int]] = []
        map_multi: Dict[Any, int] = {}

        stack = [self.root]
        while stack:
            node = stack.pop()
            d = self.dist(query_seq, node.key)

            is_candidate = (
                (not exclude_same_id) or (query_id is None) or (node.id != query_id)
            )

            if is_candidate:
                if scope == "min":
                    # Enforce radius if provided (do NOT accept candidates beyond radius)
                    if (max_radius is None) or (d <= max_radius):
                        if best_d is None or d < best_d:
                            best_d = d
                            best_ids = [node.id]
                        elif d == best_d:
                            best_ids.append(node.id)

                elif scope == "radius":
                    assert max_radius is not None
                    if d <= max_radius:
                        # collect per-neighbor distance
                        ids_multi.append(node.id)
                        dists_multi.append(int(d))
                        pairs_multi.append((node.id, int(d)))
                        map_multi[node.id] = int(d)  # dedup by id

                else:  # scope == "exact"
                    assert max_radius is not None
                    if d == max_radius:
                        ids_multi.append(node.id)
                        dists_multi.append(int(d))
                        pairs_multi.append((node.id, int(d)))
                        map_multi[node.id] = int(d)

            # --- BK-tree pruning band ---
            if max_radius is not None:
                # Tight band is fine for all scopes; candidate acceptance is filtered above
                lo = max(0, d - max_radius)
                hi = d + max_radius
            else:
                # No explicit radius: use current best distance as dynamic radius (classic BK-tree)
                radius = best_d if best_d is not None else d
                lo = max(0, d - radius)
                hi = d + radius

            for edge_d, child in node.children.items():
                if lo <= edge_d <= hi:
                    stack.append(child)

        # ---- finalize output ----
        if scope == "min":
            if not best_ids:
                # No candidate within radius (if set), or empty tree
                return ([], None) if tie == "all" else (None, None)
            # deterministic order, dedup
            best_ids = sorted(set(best_ids))
            if tie == "all":
                # For scope=min, returns "parallel"/"pairs"/"map" don't add value; keep (list_ids, dist)
                return (best_ids, best_d)
            elif tie == "min_id":
                return (min(best_ids), best_d)
            else:  # "first" -> first in sorted order
                return (best_ids[0], best_d)

        # scope in {"radius","exact"}
        if not ids_multi:
            # Empty collection for multi-hit scopes
            if returns == "parallel":
                return ([], [])
            elif returns == "pairs":
                return []
            elif returns == "map":
                return {}
            else:  # "single" makes no sense here; return empty parallel by default
                return ([], [])

        # Deduplicate + sort deterministically
        # Keep (id,dist) alignment; easiest path: unique by id using a dict, then sort by (dist, id).
        uniq = {}
        for i, d in zip(ids_multi, dists_multi):
            uniq[i] = d
        items = sorted(uniq.items(), key=lambda kv: (kv[1], str(kv[0])))

        if returns == "parallel":
            ids_sorted = [k for k, _ in items]
            dist_sorted = [int(v) for _, v in items]
            return ids_sorted, dist_sorted
        elif returns == "pairs":
            return [(k, int(v)) for k, v in items]
        elif returns == "map":
            return {k: int(v) for k, v in items}
        else:
            # returns="single" doesn't apply to multi-hit scopes -> fallback: best (min dist, then min id)
            best_k, best_v = items[0]
            return best_k, int(best_v)


# -------------------------------
# Public helpers
# -------------------------------


def build_ref_index(
    ref_ids: List[Any],
    ref_seqs: List[str],
    distance: Union[str, Callable[[str, str], Any], DistanceSpec] = "Levenshtein",
    *,
    verify: bool = False,
    verify_samples: int = 300,
    rng_seed: int = 42,
) -> BKTree:
    """
    Build and return a BKTree index over the reference set.
    Use this once, then pass the resulting 'tree' to query().

    distance:
      - string: any rapidfuzz distance name with a .distance integer metric
                (e.g., 'Levenshtein','Hamming','Indel','OSA','LCSseq', ...)
      - callable: custom (a,b)->int
      Aliases: 'lev','edit','ham'

    verify:
      - If TRUE, samples the given sequences to check integer outputs, identity,
        symmetry, and triangle inequality before building (fast sanity check).
    """
    if len(ref_ids) != len(ref_seqs):
        raise ValueError("ref_ids and ref_seqs must have the same length")

    spec = distance if isinstance(distance, DistanceSpec) else DistanceSpec(distance)

    if verify:
        spec.assert_ok(
            ref_seqs,
            n_pairs=min(verify_samples, max(50, len(ref_seqs))),
            n_triples=verify_samples,
            rng_seed=rng_seed,
        )

    tree = BKTree(distance=spec)
    tree.build(ref_seqs, ref_ids)
    return tree


def query(
    tree: BKTree,
    query_ids: Union[List[Any], Any],
    query_seqs: Union[List[str], str],
    *,
    mode: str = "single",  # "one" | "single" | "all"
    max_radius: Optional[int] = None,
    exclude_same_id: bool = True,
    scope: str = "min",  # "min" | "radius" | "exact"
    tie: str = "min_id",  # for scope="min": "first" | "min_id" | "all"
    returns: str = "parallel",  # for multi-hit scopes: "parallel" | "pairs" | "map" | "edges"
) -> QueryResult:
    """
    Unified wrapper for batch queries.
      - mode="one":    query a single (id, seq); returns one record as from tree.nearest(...)
      - mode="single": iterate over vectors; returns lists aligned with query_ids
      - mode="all":    like "single", but also supports returns="edges" to emit a flat edge list

    Notes:
      * For scope="min": 'tie' controls collapsing of ties at the minimal distance.
      * For scope in {"radius","exact"}:
          - 'returns' controls shape: "parallel", "pairs", "map", or "edges" (edges only for mode="all")
          - 'tie' is ignored (there is no "best" when you want all within the scope).
    """
    mode = mode.lower()
    scope = scope.lower()
    tie = tie.lower()
    returns = returns.lower()

    # ---- mode="one": accept scalar inputs (or length-1 lists) ----
    if mode == "one":
        if isinstance(query_ids, list):
            if len(query_ids) != 1:
                raise ValueError('mode="one" expects a single query_id')
            qid = query_ids[0]
        else:
            qid = query_ids
        if isinstance(query_seqs, list):
            if len(query_seqs) != 1:
                raise ValueError('mode="one" expects a single query_seq')
            qseq = query_seqs[0]
        else:
            qseq = query_seqs

        # For scope="min": you can use tie="all"/"first"/"min_id"; returns is ignored.
        # For multi-hit scopes: returns controls shape; tie ignored.
        if scope == "min":
            return cast(
                MinOneResult | MinAllResult,
                tree.nearest(
                    qseq,
                    qid,
                    max_radius=max_radius,
                    exclude_same_id=exclude_same_id,
                    tie=tie,
                    scope=scope,
                    returns="single",
                ),
            )
        raw = tree.nearest(
            qseq,
            qid,
            max_radius=max_radius,
            exclude_same_id=exclude_same_id,
            tie="all",
            scope=scope,
            returns=returns,
        )
        if returns == "parallel":
            return cast(OneParallelResult, raw)
        if returns == "pairs":
            return cast(OnePairsResult, raw)
        if returns == "map":
            return cast(OneMapResult, raw)
        return cast(QueryResult, raw)

    # ---- modes requiring vector inputs ----
    if not (isinstance(query_ids, list) and isinstance(query_seqs, list)):
        raise ValueError(f'mode="{mode}" expects query_ids and query_seqs as lists')
    if len(query_ids) != len(query_seqs):
        raise ValueError("query_ids and query_seqs must have the same length")

    n = len(query_ids)

    # scope="min": 'single' means one best per query (id + dist)
    if scope == "min" and mode == "single":
        nn_ids: List[Optional[Any]] = [None] * n
        nn_dsts: List[Optional[int]] = [None] * n
        for i, (qid, qseq) in enumerate(zip(query_ids, query_seqs)):
            ans_id, ans_d = cast(
                MinOneResult,
                tree.nearest(
                    qseq,
                    qid,
                    max_radius=max_radius,
                    exclude_same_id=exclude_same_id,
                    tie=tie,
                    scope=scope,
                    returns="single",
                ),
            )
            nn_ids[i] = ans_id
            nn_dsts[i] = ans_d if ans_d is not None else None
        return nn_ids, nn_dsts

    # For multi-hit scopes ("radius"/"exact") OR mode="all":
    # returns = "parallel" | "pairs" | "map" | "edges"
    if returns == "edges":
        # Build a flat edge list: (query_id, neighbor_id, distance) across all queries
        edges: List[Tuple[Any, Any, int]] = []
        for qid, qseq in zip(query_ids, query_seqs):
            ids_i, dists_i = cast(
                OneParallelResult,
                tree.nearest(
                    qseq,
                    qid,
                    max_radius=max_radius,
                    exclude_same_id=exclude_same_id,
                    tie="all",
                    scope=scope,
                    returns="parallel",
                ),
            )
            for nid, dist in zip(ids_i, dists_i):
                edges.append((qid, nid, int(dist)))
        return edges

    # Otherwise return per-query collections (list-columns friendly)
    if returns == "parallel":
        all_ids: List[List[Any]] = [[] for _ in range(n)]
        all_dist: List[List[int]] = [[] for _ in range(n)]
        for i, (qid, qseq) in enumerate(zip(query_ids, query_seqs)):
            ids_i, dists_i = cast(
                OneParallelResult,
                tree.nearest(
                    qseq,
                    qid,
                    max_radius=max_radius,
                    exclude_same_id=exclude_same_id,
                    tie="all",
                    scope=scope,
                    returns="parallel",
                ),
            )
            all_ids[i] = ids_i
            all_dist[i] = dists_i
        return all_ids, all_dist

    elif returns == "pairs":
        all_pairs: List[List[Tuple[Any, int]]] = [[] for _ in range(n)]
        for i, (qid, qseq) in enumerate(zip(query_ids, query_seqs)):
            pairs_i = cast(
                OnePairsResult,
                tree.nearest(
                    qseq,
                    qid,
                    max_radius=max_radius,
                    exclude_same_id=exclude_same_id,
                    tie="all",
                    scope=scope,
                    returns="pairs",
                ),
            )
            all_pairs[i] = pairs_i
        return all_pairs

    elif returns == "map":
        all_maps: List[Dict[Any, int]] = [{} for _ in range(n)]
        for i, (qid, qseq) in enumerate(zip(query_ids, query_seqs)):
            map_i = cast(
                OneMapResult,
                tree.nearest(
                    qseq,
                    qid,
                    max_radius=max_radius,
                    exclude_same_id=exclude_same_id,
                    tie="all",
                    scope=scope,
                    returns="map",
                ),
            )
            all_maps[i] = map_i
        return all_maps

    else:
        raise ValueError(
            'returns must be one of "parallel", "pairs", "map", or "edges"'
        )


# %% Usage Notes

"""
This module has no built-in argparse interface yet; the snippets above show
idiomatic "CLI" usage via `python -c` or here-doc scripts.

USAGE (from a terminal / Python CLI)
====================================

This module exposes a BK-tree index and helpers for nearest-neighbor search
under edit-like distances (RapidFuzz backend). It is designed to be imported
and called from the Python CLI (e.g., `python -c "..."`) or tiny scripts.

Quick import check
------------------
$ python -c "import bk_tree as b; print('OK, versionless module loaded')"

Build an index and query one sequence (minimal distance)
--------------------------------------------------------
$ python - <<'PY'
from bk_tree import build_ref_index, query
ref_ids  = ['r1','r2','r3','r4']
ref_seqs = ['apple','apply','maple','snapple']
tree = build_ref_index(ref_ids, ref_seqs, distance='Levenshtein', verify=False)
(qid, qseq) = ('q1','appl')
ans_id, ans_d = query(tree, [qid], [qseq], mode='one', scope='min', tie='min_id')
print({'query_id': qid, 'best_id': ans_id, 'dist': ans_d})
PY

Batch: one best neighbor per query (vectorized)
-----------------------------------------------
$ python - <<'PY'
from bk_tree import build_ref_index, query
ref_ids  = [f'r{i}' for i in range(5)]
ref_seqs = ['cat','cot','coat','scat','at']
tree = build_ref_index(ref_ids, ref_seqs, distance='Levenshtein')
q_ids  = ['q1','q2','q3']
q_seqs = ['cut','coats','cats']
nn_ids, nn_d = query(tree, q_ids, q_seqs, mode='single', scope='min', tie='min_id')
for i,(qi,qs) in enumerate(zip(q_ids,q_seqs)):
    print(qi, qs, '->', nn_ids[i], nn_d[i])
PY

Radius search: all neighbors within a distance
----------------------------------------------
$ python - <<'PY'
from bk_tree import build_ref_index, query
ref_ids  = ['a','b','c','d']
ref_seqs = ['GATTACA','GACTATA','GATTTCA','CAT']
tree = build_ref_index(ref_ids, ref_seqs)
q_ids  = ['q1','q2']
q_seqs = ['GATTATA','GATACA']
all_ids, all_d = query(tree, q_ids, q_seqs,
                       mode='single', scope='radius',
                       max_radius=2, returns='parallel')
for qi, ids, ds in zip(q_ids, all_ids, all_d):
    print(qi, list(zip(ids, ds)))
PY

Exact-distance matches only
---------------------------
$ python - <<'PY'
from bk_tree import build_ref_index, query
tree = build_ref_index(['r1','r2'], ['foo','foe'])
print(query(tree, ['q1'], ['fob'], mode='single',
            scope='exact', max_radius=1, returns='pairs'))
PY

Edge list for graph workflows (across all queries)
--------------------------------------------------
$ python - <<'PY'
from bk_tree import build_ref_index, query
tree = build_ref_index(['r1','r2','r3'], ['AAA','AAT','TAA'])
edges = query(tree,
              ['q1','q2'], ['AAT','TAT'],
              mode='all', scope='radius',
              max_radius=1, returns='edges')
# Emits (query_id, neighbor_id, distance) tuples
for e in edges: print(e)
PY

Using alternative distances (Hamming, etc.)
-------------------------------------------
# Hamming requires equal-length strings:
$ python - <<'PY'
from bk_tree import build_ref_index, query
tree = build_ref_index(['r1','r2','r3'], ['ACGT','ACGA','TCGT'], distance='Hamming')
print(query(tree, ['q1'], ['ACGG'], mode='one', scope='min'))
PY

Light sanity check of the distance metric (optional)
----------------------------------------------------
$ python - <<'PY'
from bk_tree import build_ref_index
# verify=True samples the sequences to check identity/symmetry/triangle inequality
tree = build_ref_index(['r1','r2','r3'], ['apple','apply','maple'],
                       distance='Levenshtein', verify=True, verify_samples=200)
print('verified index ready')
PY

Minimal R/reticulate example!
-----------------------------
# Need this
library(reticulate)

# Point path to where the py script lives
bk <- import_from_path("bk_tree", path = "./py/")

ids <- c("SIVmac239M2.111392","SIVmac239M2.630","SIVmac239M2.83869")                                                      #a list of barcode names
seqs <- c("ACGCGCGCATCGGTGCGTATCGCATGGCTCGCGT","ACGCGCGCATCGATTTTTACTGCATGGCTCGCGT","ACGCGCGCATCGGCAGGGAAAACATGGCTCGCGT") #a list of barcode seqs in the same order

# Build your reference bk_tree
ref_tree <- bk$build_ref_index(ids, seqs)

# Query the tree
res <- bk$query(
  tree = ref_tree,
  query_ids = ids,
  query_seqs = seqs,
  mode = "all",
  tie = "all",
  scope = "radius",
  max_radius = 2,
  exclude_same_id = T,
  returns = "edges"
)
"""

#%% Version history

"""
1.0.1 - Minor docstring and type annotation improvements.
1.0.0 - Initial release of bk_tree.py with BK-tree implementation, distance metric support, and query interface.
"""
