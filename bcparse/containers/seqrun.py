# seqrun.py

"""
Name:       seqrun.py
Author:     CAG
Version:    1.1.0
Date:       2026/03/26
Refactored: 2026/05/26

`SeqRun` is the main in-memory representation for one sequencing run:
- many `SeqSamp` members
- normalized run-level metadata
- optional RunInfo provenance
"""

# %% Imports

from __future__ import annotations

from pathlib import Path

import pandas as pd

from bcparse.containers.runinfo import RunInfo
from bcparse.containers.seqsamp import SeqSamp

# %% SeqRun model


class SeqRun:
    """
    Strict container for one sequencing run made up of SeqSamp objects.

    Construction is owned by parse/compile mode code. SeqRun only validates and
    stores the already-materialized samples, run metadata, optional RunInfo, and
    canonical long dataframe.
    """

    def __init__(
        self,
        samples: dict[str, SeqSamp],
        *,
        run_number,
        run_name,
        run_date,
        run_id,
        runinfo: RunInfo | None,
        filepath,
        source="unknown",
        data_df: pd.DataFrame,
    ):
        self._validate_samples(samples)
        if data_df.empty:
            raise ValueError("SeqRun requires a non-empty data_df.")

        self.samples = dict(samples)
        self.run_number = run_number
        self.run_name = run_name
        self.run_date = run_date
        self.run_id = run_id
        self.runinfo = runinfo
        self.filepath = Path(filepath) if filepath else None
        self.source = source
        self.data_df = data_df.copy()

    @staticmethod
    def _validate_samples(samples: dict[str, SeqSamp]) -> None:
        if not samples:
            raise ValueError("SeqRun requires at least one SeqSamp.")
        for sample_id, sample in samples.items():
            resolved_id = str(sample_id)
            if sample.sample_id != resolved_id:
                raise ValueError(
                    f"SeqRun key/sample_id mismatch: key={resolved_id}, object={sample.sample_id}"
                )

    def ensure_runinfo(self, filepath: str | Path | None = None) -> RunInfo:
        if self.runinfo is not None:
            return self.runinfo
        self.runinfo = RunInfo.from_long_df(self.df, filepath=filepath)
        return self.runinfo

    # ====================
    # EMIT
    # ====================

    @property
    def sample_list(self) -> list[SeqSamp]:
        return list(self.samples.values())

    @property
    def sample_ids(self) -> list[str]:
        return list(self.samples.keys())

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    @property
    def sample_meta_df(self) -> pd.DataFrame:
        return pd.DataFrame([s.sample_meta for s in self.sample_list])

    @property
    def df(self) -> pd.DataFrame:
        return self.data_df.copy()

    @property
    def df_above_cutoff(self) -> pd.DataFrame:
        df = self.df
        if "above_cutoff" not in df.columns:
            return df.iloc[0:0].copy()
        return df.loc[df["above_cutoff"].fillna(False)].copy()

    @property
    def groups(self) -> list[str]:
        df = self.df
        if "samp_group" not in df.columns:
            return []
        return df["samp_group"].dropna().astype(str).unique().tolist()


# %% Usage notes
"""
SeqRun is the per-run container used by both parse and compile modes. It
represents one sequencing run after ingest has already produced SeqSamp
children and a run-scoped normalized long-format dataframe.

Construction:
 - ingest/parse_mode.py builds one SeqRun directly at the end of ParseMode.
 - ingest/compile_mode.py builds one SeqRun per run_id while materializing the
   final RunSeries.
 - SeqRun does not own dataframe-to-run construction helpers; mode classes
   decide how samples are grouped and how run metadata is recovered.

Stored state:
 - samples is the canonical {sample_id: SeqSamp} mapping.
 - data_df is the cached run-scoped long dataframe used by the df property.
 - run_number, run_name, run_date, and run_id identify the sequencing run.
 - runinfo stores workbook-derived or synthesized RunInfo metadata.
 - filepath and source record ingest provenance.

Contracts:
 - one SeqRun maps to exactly one run_id.
 - samples must be non-empty and keyed by each SeqSamp.sample_id.
 - data_df must be non-empty.
 - runinfo may be None during compile-mode construction for base-CSV-derived
   runs, then filled by ensure_runinfo().

Downstream consumers:
 - RunSeries stores SeqRun objects in its runs mapping.
 - emit/workbooks.py uses df, sample_list, groups, and runinfo for parse and
   compile workbook output.
 - emit/views.py uses sample_list and df for grouped output matrices and
   concatenated runinfo views.
"""


# %% Versions
"""
v1.1.0 20260625
- Removed from_long_df() and flexible sample coercion.
- ParseMode and CompileMode now own mode-specific SeqRun materialization.

v1.0.0 20260526
- Initial version, accomodating new model/container structure
- Supercedes deprecated parse_countdata module and associated code in bcParse.py
"""
