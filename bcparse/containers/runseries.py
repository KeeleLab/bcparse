# runseries.py
"""
Name:       runseries.py
Author:     CAG
Version:    2.1.0
Date:       2026/06/23

`RunSeries` is the top-level model used by compile mode, holding a series of `SeqRun` objects.
"""

# %% Imports

from __future__ import annotations

from typing import Any

import pandas as pd

from bcparse.containers.runinfo import RunInfo
from bcparse.containers.seqrun import SeqRun

# %% RunSeries model


class RunSeries:
    """
    Container for a series of sequencing runs.
    Holds multi-run compiled data and per-run metadata for compile mode.

    Identity contracts for compile mode:
    - samples are identified by `sample_id`; composite group::name::date::run_number
    - runs are identified by `run_id`; composite run_name::run_number
    """

    def __init__(
        self,
        runs: dict[str, SeqRun],
        series_meta: dict[str, Any],
        data_df: pd.DataFrame,
        runinfo_by_run_id: dict[str, RunInfo],
        source: str = "unknown",
    ):
        self.runs = {str(run_id): run for run_id, run in runs.items()}
        self.series_meta = dict(series_meta)
        self.data_df = data_df.copy()
        self.runinfo_by_run_id = {str(k): v for k, v in runinfo_by_run_id.items()}
        self.source = source

    # ====================
    # INGEST
    # ====================

    @classmethod
    def from_long_df(
        cls,
        df: pd.DataFrame,
        *,
        runinfo_by_run_id: dict[str, RunInfo],
        source: str = "unknown",
    ) -> "RunSeries":
        """
        Materialize a RunSeries from a normalized long-format dataframe.
        Splits by run identity and builds one SeqRun per run.
        """
        ri = runinfo_by_run_id
        # Base CSVs can pass an empty runinfo map; compile.py fills it later
        # with ensure_runinfo().

        if df.empty:  # safeguard if called directly for some reason
            raise ValueError("RunSeries.from_long_df requires a non-empty dataframe.")

        # Split by run_id and build SeqRun objects
        runs = {
            str(run_id): SeqRun.from_long_df(
                run_df,
                source=source,
                runinfo=ri.get(str(run_id)),
            )
            for run_id, run_df in df.groupby("run_id", sort=False)
        }

        # Return a RunSeries object with the runs, series metadata, and the original dataframe
        return cls(
            runs=runs,
            series_meta={"groupby_col": "run_id"},
            data_df=df,
            runinfo_by_run_id=ri,
            source=source,
        )

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
        return self.data_df.copy()


# %% Repository usage notes
"""
RunSeries is the compile-mode top-level container. It represents a
deduplicated series of sequencing runs after compile.py has normalized,
merged, and QC-annotated long-format input data.

Construction:
 - ingest/compile.py is the only current external construction path.
 - build_runseries() creates a final combined long dataframe, then calls
   RunSeries.from_long_df().
 - from_long_df() splits the dataframe by run_id and builds one SeqRun per run.
 - Empty dataframes are invalid here; compile.py should fail before calling
   this constructor if no usable input data exists.

Stored state:
 - runs is the canonical {run_id: SeqRun} mapping.
 - data_df is the cached compiled long dataframe used by the df property.
 - series_meta stores compile-level metadata such as groupby_col, out_prefix,
   and n_source_runs.
 - runinfo_by_run_id stores per-run RunInfo objects after compile finalization.

Downstream consumers:
 - emit/workbooks.py uses RunSeries for compile workbook and CSV output.
 - emit/workbooks.py::emit_qc() consumes run_list or df depending on whether
   collapse_to_parent is enabled.
 - emit/views.py::concat_series_runinfo() walks run_list and concatenates each
   run's raw runinfo table for compile output.
"""

# %% Versions
"""
v2.1.0 20260623
 - Dropped over-flexible input types in lieu of trusting upstream contract
 - Formatting and comments

v2.0.0 20260526
 - Plain class replacing dataclass
 - copy() and collapse_samples_in_place() removed; emit_qc works directly from sample_list
 - Compile path collapsed from ~15 staticmethods to build_from_parsed_files + _dedup + _core_qc
 - sample_id is the compile deduplication contract
 - Refactored into bcparse package structure

v1.0.0
 - Initial version, accomodating new model/container structure
 - Supercedes deprecated xlsx_compiler flow pre-v4
"""
