# seqsamp.py
"""
Name:       seqsamp.py
Author:     CAG
Version:    1.0
Date:       2026/03/26

Per-sample data model.

`SeqSamp` is the smallest biological analysis unit in the refactor:
- one logical sample
- one barcode/count dataframe scoped to that sample
- one metadata dictionary describing that sample
"""

#%% Imports

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from utils.runinfo import make_sample_id

#%% SeqSamp model

@dataclass
class SeqSamp:
    """
    Canonical per-sample contract for downstream analysis.

    Important invariant: one `SeqSamp` must map to exactly one logical sample.
    Many later operations, especially collapse/QC, rely on that boundary.

    In normal parse/workbook-backed flows that sample identity is `idx_name`.
    Composite ids are retained only as a fallback for rows that do not carry an
    idx_name, such as some legacy/base-csv inputs.
    """

    
    df:          pd.DataFrame                                 # Long-format barcode table restricted to one logical sample.
    sample_meta: dict[str, Any] = field(default_factory=dict) # Metadata for the same logical sample; copied onto `df` when missing.
    sample_id:   str | None = None                            # Stable sample identifier shared across parse and compile paths.
    source:      str = "unknown"                              # Provenance tag for debugging / reconstruction.

    def __post_init__(self) -> None:
        # Work on private copies so model instances are safe to reuse.
        self.df = self.df.copy()
        self.sample_meta = dict(self.sample_meta)
        df_first = self.df.iloc[0] if not self.df.empty else None
        row_sample_ids: list[str] = []

        # If the dataframe already carries sample ids, keep them available as a
        # preferred source of truth for compiled/base-data round-trips.
        if "sample_id" in self.df.columns:
            row_sample_ids = [
                value
                for value in self.df["sample_id"].dropna().astype(str).unique().tolist()
                if value
            ]

        # Pull obvious sample-level fields off the dataframe when explicit
        # metadata was not provided by the caller.
        if df_first is not None:
            for key in ["idx_name", "samp_group", "samp_name", "samp_date", "input", "run_number"]:
                if key not in self.sample_meta and key in self.df.columns:
                    value = df_first.get(key)
                    if pd.notna(value):
                        self.sample_meta[key] = value

        # Resolve `sample_id` from strongest to weakest available evidence.
        # Prefer idx_name whenever it exists; that is the intended sample
        # identity within a run.
        if self.sample_id is None:
            self.sample_id = (
                self.sample_meta.get("sample_id")
                or (row_sample_ids[0] if len(row_sample_ids) == 1 else None)
                or self.sample_meta.get("idx_name")
                or (df_first.get("idx_name") if df_first is not None else None)
                or make_sample_id(
                    self.sample_meta.get("samp_group"),
                    self.sample_meta.get("samp_name"),
                    self.sample_meta.get("samp_date"),
                    self.sample_meta.get("run_number"),
                )
                or (
                    make_sample_id(
                        df_first.get("samp_group"),
                        df_first.get("samp_name"),
                        df_first.get("samp_date"),
                        df_first.get("run_number"),
                    )
                    if df_first is not None
                    else None
                )
                or self.sample_meta.get("samp_name")
                or (df_first.get("samp_name") if df_first is not None else None)
            )

        if self.sample_id is None:
            raise ValueError("SeqSamp requires a resolvable sample_id.")

        self.sample_id = str(self.sample_id)
        self.sample_meta["sample_id"] = self.sample_id

        # Ensure dataframe and metadata stay aligned on the same sample fields.
        for key, value in self.sample_meta.items():
            if key not in self.df.columns:
                self.df[key] = value

        if "sample_id" not in self.df.columns:
            self.df["sample_id"] = self.sample_id

        # Hard-stop if this object accidentally spans more than one sample.
        self._validate_single_sample_scope()

    def _validate_single_sample_scope(self) -> None:
        """
        Enforce the core invariant that one SeqSamp maps to one logical sample.

        Note that rows may still imply multiple metadata-derived composite ids
        in legacy data, but when idx_name is present it remains the operative
        identity for this object.
        """
        if self.df.empty:
            return

        sample_ids = self.df["sample_id"].dropna().astype(str).unique().tolist()
        if len(sample_ids) > 1:
            raise ValueError(
                f"SeqSamp expected one sample_id, found {sample_ids}"
            )
        if sample_ids and sample_ids[0] != self.sample_id:
            raise ValueError(
                f"SeqSamp sample_id mismatch: object={self.sample_id}, rows={sample_ids[0]}"
            )

        run_col = "run_number"
        req_cols = {"samp_group", "samp_name", "samp_date", run_col}
        if not req_cols.issubset(self.df.columns):
            return

        inferred_ids = (
            self.df.apply(
                lambda row: make_sample_id(
                    row.get("samp_group"),
                    row.get("samp_name"),
                    row.get("samp_date"),
                    row.get(run_col),
                ),
                axis=1,
            )
            .dropna()
            .astype(str)
        )

        inferred_ids = [value for value in inferred_ids.unique().tolist() if value != ":::"]
        if len(inferred_ids) > 1:
            raise ValueError(
                f"SeqSamp rows span multiple logical samples: {inferred_ids}"
            )

    @classmethod
    def from_table(
        cls,
        *,
        sample_df: pd.DataFrame,
        sample_meta: dict[str, Any] | None = None,
        sample_id: str | None = None,
        source: str = "unknown",
    ) -> "SeqSamp":
        """
        Convenience constructor for turning one sample-scoped dataframe into a
        model, preferring idx_name as the resolved sample identity.
        """
        meta = {} if sample_meta is None else dict(sample_meta)
        first_row = sample_df.iloc[0] if not sample_df.empty else None

        # Rehydrate metadata from the table itself when round-tripping from CSV /
        # combined long data.
        if first_row is not None:
            for key in ["sample_id", "idx_name", "samp_group", "samp_name", "samp_date", "input", "run_number"]:
                if key not in meta and key in sample_df.columns:
                    value = first_row.get(key)
                    if pd.notna(value):
                        meta[key] = value

        if sample_id is None:
            sample_id = (
                meta.get("sample_id")
                or meta.get("idx_name")
                or make_sample_id(
                    meta.get("samp_group"),
                    meta.get("samp_name"),
                    meta.get("samp_date"),
                    meta.get("run_number"),
                )
                or "::".join(
                    str(meta.get(col, ""))
                    for col in ["run_number", "filename", "samp_name", "samp_date"]
                    if col in meta
                )
                or "::".join(
                    str(meta.get(col, ""))
                    for col in ["samp_group", "samp_name", "samp_date"]
                    if col in meta
                )
            )

        if sample_id is not None:
            meta["sample_id"] = sample_id

        return cls(
            df=sample_df.reset_index(drop=True),
            sample_meta=meta,
            sample_id=sample_id,
            source=source,
        )

    def collapse_to_parent(self) -> "SeqSamp":
        """
        Return a copy of this SeqSamp with child rows collapsed to parent rows.

        This is the safe/default collapse scope for view generation because it
        constrains qc.collapse_to_parent() to one sample at a time.
        """
        import utils.qc_ops as qc

        # Collapse within the sample boundary, then rebuild another SeqSamp so
        # downstream code keeps the same object contract.
        collapsed_df = qc.collapse_to_parent(self.df.copy())
        if isinstance(collapsed_df, tuple):
            collapsed_df = collapsed_df[0]

        return SeqSamp(
            df=collapsed_df,
            sample_meta=self.sample_meta.copy(),
            sample_id=self.sample_id,
            source=self.source,
        )

#%% Versions
"""
v1.0 
- Initial version, accomodating new model/container structure
- Support structure for seqrun and runseries, this is the planned new atomic unit for downstream analysis.
"""
