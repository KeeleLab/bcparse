# seqsamp.py
"""
Name:      seqsamp.py
Author:     CAG
Version:    1.0
Date:       2026/03/26
Refactored: 2026/05/26

Per-sample data model.

`SeqSamp` is the smallest biological analysis unit in the refactor:
- one logical sample
- one barcode/count dataframe scoped to that sample
- one metadata dictionary describing that sample
"""

#%% Imports

from __future__ import annotations

from typing import Any

import pandas as pd

#%% SeqSamp model

class SeqSamp:
    """
    Canonical per-sample contract for downstream analysis.

    Important invariant: one `SeqSamp` must map to exactly one logical sample.
    Many later operations, especially collapse/QC, rely on that boundary.

    The caller must provide sample_id directly or in sample_meta.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        sample_meta: dict[str, Any] | None = None,
        sample_id: str | None = None,
        source: str = "unknown",
    ):
        # Work on private copies so model instances are safe to reuse.
        self.df = df.copy()
        self.sample_meta = dict(sample_meta or {})
        self.sample_id = sample_id
        self.source = source

        if self.sample_id is None:
            self.sample_id = self.sample_meta.get("sample_id")

        if self.sample_id is None or pd.isna(self.sample_id) or str(self.sample_id).strip() == "":
            raise ValueError("SeqSamp requires sample_id.")

        self.sample_id = str(self.sample_id)
        self.sample_meta["sample_id"] = self.sample_id

        # Ensure dataframe and metadata stay aligned on the same sample fields.
        for key, value in self.sample_meta.items():
            if key not in self.df.columns:
                self.df[key] = value

        if "sample_id" not in self.df.columns:
            self.df["sample_id"] = self.sample_id

        self._validate_single_sample_scope()

    def _validate_single_sample_scope(self) -> None:
        """
        Enforce the core invariant that one SeqSamp maps to one logical sample.

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
        Convenience constructor for turning one sample-scoped dataframe into a model.
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
            sample_id = meta.get("sample_id")

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
        from bcparse.qc.collapse import collapse_to_parent

        # Collapse within the sample boundary, then rebuild another SeqSamp so
        # downstream code keeps the same object contract.
        collapsed_df = collapse_to_parent(self.df.copy())
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
2026-05-26 - Refactored into bcparse package structure
"""
