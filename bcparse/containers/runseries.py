# runseries.py
"""
Name:       runseries.py
Author:     CAG
Version:    2.0
Date:       2026/05/21
Refactored: 2026/05/26

`RunSeries` is the top-level model used by compile mode:
- n `SeqRun` objects + a cached raw compiled long dataframe
- round-trip construction from canonical long-format data
"""

#%% Imports

from __future__ import annotations

from typing import Any, Iterable, Mapping

import pandas as pd

from bcparse.containers.runinfo import RunInfo
from bcparse.containers.seqrun import SeqRun

#%% RunSeries model

class RunSeries:
    """
    Container for a series of sequencing runs.
    Holds multi-run compiled data and per-run metadata for compile mode.

    Identity contract for compile mode:
    - samples are identified by `sample_id`
    - `sample_id` is the composite group/name/date/run_number identity
    - `idx_name` remains index metadata, not the sample identity contract
    """

    def __init__(
        self,
        runs: dict[str, SeqRun] | Iterable[SeqRun] | None = None,
        series_meta: dict[str, Any] | None = None,
        data_df: pd.DataFrame | None = None,
        runinfo_by_run_id: dict[str, RunInfo] | None = None,
        source: str = "unknown",
    ):
        self.runs             = self._coerce_runs(runs or {})
        self.series_meta      = dict(series_meta or {})
        self.data_df          = None if data_df is None else data_df.copy()
        self.runinfo_by_run_id = {
            str(k): v for k, v in (runinfo_by_run_id or {}).items() if v is not None
        }
        self.source = source

    # ====================
    # INGEST
    # ====================

    @classmethod
    def from_long_df(
        cls,
        df: pd.DataFrame,
        *,
        source: str = "unknown",
        runinfo_by_run_id: dict[str, RunInfo] | None = None,
    ) -> "RunSeries":
        """
        Materialize a RunSeries from a normalized long-format dataframe.
        Splits by run identity and builds one SeqRun per run.
        """
        ri = runinfo_by_run_id or {}

        if df.empty:
            return cls(runinfo_by_run_id=ri, source=source)

        runs = {
            str(run_id): SeqRun.from_long_df(
                run_df,
                source=source,
                runinfo=ri.get(str(run_id)),
            )
            for run_id, run_df in df.groupby("run_id", sort=False)
        }

        return cls(
            runs=runs,
            series_meta={"groupby_col": "run_id"},
            runinfo_by_run_id=ri,
            source=source,
        )

    @staticmethod
    def _coerce_runs(runs: Mapping[str, SeqRun] | Iterable[SeqRun]) -> dict[str, SeqRun]:
        """Normalize runs to {run_id: SeqRun}."""
        items = (
            list(runs.items())
            if isinstance(runs, Mapping)
            else [(run.run_id, run) for run in runs]
        )
        run_ids = [str(run_id) for run_id, _ in items]

        dup_ids = sorted({run_id for run_id in run_ids if run_ids.count(run_id) > 1})
        if dup_ids:
            raise ValueError(f"Duplicate run_id in RunSeries: {', '.join(dup_ids)}")

        return dict(zip(run_ids, [run for _, run in items]))

    # ====================
    # EMIT
    # ====================

    @property
    def run_list(self) -> list[SeqRun]:
        return list(self.runs.values())

    @property
    def run_ids(self) -> list[str]:
        return list(self.runs.keys())

    @property
    def n_runs(self) -> int:
        return len(self.runs)

    @property
    def df(self) -> pd.DataFrame:
        if self.data_df is not None:
            return self.data_df.copy()
        if not self.runs:
            return pd.DataFrame()
        return pd.concat([run.df for run in self.run_list], ignore_index=True)

#%% Versions
"""
v1.0 
- Initial version, accomodating new model/container structure
- Supercedes deprecated xlsx_compiler module and associated code in bcParse.py

v2.0
- Plain class replacing dataclass
- copy() and collapse_samples_in_place() removed; emit_qc works directly from sample_list
- Compile path collapsed from ~15 staticmethods to build_from_parsed_files + _dedup + _core_qc
- sample_id is the compile deduplication contract
2026-05-26 - Refactored into bcparse package structure
"""
