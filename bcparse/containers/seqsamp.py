# seqsamp.py
"""
Name:      seqsamp.py
Author:     CAG
Version:    1.1.0
Date:       20260625

Per-sample data model.

`SeqSamp` is the smallest biological analysis unit in the refactor:
- one logical sample
- one barcode/count dataframe scoped to that sample
- one metadata dictionary describing that sample
"""

# %% Imports

from __future__ import annotations

from typing import Any, cast

import pandas as pd

# %% SeqSamp model


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

        if (
            self.sample_id is None
            or pd.isna(self.sample_id)
            or str(self.sample_id).strip() == ""
        ):
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
            raise ValueError(f"SeqSamp expected one sample_id, found {sample_ids}")
        if sample_ids and sample_ids[0] != self.sample_id:
            raise ValueError(
                f"SeqSamp sample_id mismatch: object={self.sample_id}, rows={sample_ids[0]}"
            )

    def collapse_to_parent(self) -> "SeqSamp":
        """
        Return a copy of this SeqSamp with child rows collapsed to parent rows.

        This is the safe/default collapse scope for view generation because it
        constrains qc.collapse_to_parent() to one sample at a time.

        Currently, we only collapse at the sample level,
        and doing this here enforces that contract.
        """
        from bcparse.qc.collapse import collapse_to_parent

        # Collapse within the sample boundary, then rebuild another SeqSamp so
        # downstream code keeps the same object contract.
        collapsed_df = collapse_to_parent(self.df.copy())
        if isinstance(collapsed_df, tuple):
            collapsed_df = collapsed_df[0]
        collapsed_df = cast(pd.DataFrame, collapsed_df)

        return SeqSamp(
            df=collapsed_df,
            sample_meta=self.sample_meta.copy(),
            sample_id=self.sample_id,
            source=self.source,
        )


# %% Usage notes
"""
SeqSamp is the per-sample container used by both parse and compile modes.
It represents one logical sample after ingest has already produced a
sample-scoped long-format dataframe.

Construction:
 - ingest/parse_mode.py builds SeqSamp objects directly from parsed FASTQ data
   after runinfo merge and per-sample QC.
 - ingest/compile_mode.py builds SeqSamp objects directly from final compiled
   long-dataframe slices grouped by sample_id.
 - SeqSamp does not own dataframe-to-sample construction helpers; mode classes
   decide how sample metadata is recovered from their own input formats.

Stored state:
 - df is the sample-scoped long dataframe.
 - sample_meta is the metadata dictionary carried with the sample.
 - sample_id is the canonical logical sample identity.
 - source records the ingest path that built the object.

Contracts:
 - one SeqSamp maps to exactly one sample_id.
 - sample_id may be passed directly or through sample_meta, but it must resolve
   to a non-empty value.
 - constructor validation rejects dataframe rows from more than one sample_id.

Downstream consumers:
 - SeqRun stores SeqSamp objects in its samples mapping.
 - emit/workbooks.py and emit/views.py call collapse_to_parent() when output
   should collapse child barcodes within each sample.
"""


# %% Versions
"""
v1.1.0 20260625
- Tightened SeqSamp toward a light per-sample container.
- Removed from_table(); parse/compile modes now own construction from tabular
  inputs and call SeqSamp(...) directly.
- Kept _validate_single_sample_scope() as the container-level invariant check.
- Kept collapse_to_parent() as sample-scoped behavior so emit/QC collapse stays
  constrained to one logical sample.
- Applied comment/section formatting cleanup relative to the previous git
  version.

v1.0.0 20260526
- Initial version, accomodating new model/container structure
- Support structure for seqrun and runseries, this is the planned new atomic unit for downstream analysis.
"""
