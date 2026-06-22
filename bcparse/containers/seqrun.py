# seqrun.py

"""
Name:       seqrun.py
Author:     CAG
Version:    1.0
Date:       2026/03/26
Refactored: 2026/05/26

`SeqRun` is the main in-memory representation for one sequencing run:
- many `SeqSamp` members
- normalized run-level metadata
- round-trip construction from canonical long-format data
"""

#%% Imports

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

from bcparse.qc import identity as qc
from bcparse.containers.runinfo import RunInfo
from bcparse.containers.seqsamp import SeqSamp

#%% SeqRun model

class SeqRun:
    """
    Container for one sequencing run made up of arbitrary samples.
    """

    def __init__(
        self,
        samples=None,
        *,
        run_number=None,
        run_name=None,
        run_date=None,
        run_id=None,
        runinfo=None,
        filepath=None,
        source="unknown",
        data_df=None,
    ):
        self.samples    = self._coerce_samples(samples or {})
        self.run_number = run_number
        self.run_name   = run_name
        self.run_date   = run_date
        self.run_id     = run_id or qc.make_run_id(run_number, filepath, filepath)
        self.runinfo    = runinfo
        self.filepath   = Path(filepath) if filepath else None
        self.source     = source
        self.data_df    = None if data_df is None else data_df.copy()

    # ====================
    # INGEST
    # ====================

    @classmethod
    def from_long_df(
        cls,
        df: pd.DataFrame,
        *,
        source: str = "unknown",
        runinfo: RunInfo | None = None,
    ) -> "SeqRun":
        """
        Materialize a SeqRun from a normalized long-format dataframe representing one run.
        """
        if df.empty:
            return cls(source=source)

        run_ids = df["run_id"].dropna().astype(str).unique().tolist()
        if len(run_ids) > 1:
            raise ValueError(f"SeqRun.from_long_df expected one run_id, found {run_ids}")

        def first_val(col):
            return None if col not in df.columns or df[col].dropna().empty else df[col].dropna().iloc[0]

        samples = {}
        for sample_id, sample_df in df.groupby("sample_id", dropna=False, sort=False):
            sample_meta = {
                col: sample_df[col].dropna().iloc[0]
                for col in [
                    "sample_id", "idx_name", "samp_group", "samp_name", "samp_date",
                    "input", "run_number", "run_name", "run_date", "filename",
                ]
                if col in sample_df.columns and not sample_df[col].dropna().empty
            }
            resolved_id = None if pd.isna(sample_id) else str(sample_id)
            resolved_id = resolved_id or sample_meta.get("sample_id") or str(len(samples))
            samples[resolved_id] = SeqSamp.from_table(
                sample_df=sample_df,
                sample_meta=sample_meta,
                sample_id=resolved_id,
                source=source,
            )

        return cls(
            samples=samples,
            run_number=first_val("run_number"),
            run_name=first_val("run_name"),
            run_date=first_val("run_date"),
            run_id=first_val("run_id"),
            runinfo=runinfo,
            filepath=first_val("filename"),
            source=source,
        )

    @staticmethod
    def _coerce_samples(
        samples: Mapping[str, SeqSamp] | Iterable[SeqSamp],
    ) -> dict[str, SeqSamp]:
        items = samples.items() if isinstance(samples, Mapping) else ((None, s) for s in samples)
        sample_map: dict[str, SeqSamp] = {}
        for key, sample in items:
            resolved_id = str(key or sample.sample_id or sample.sample_meta.get("sample_id"))
            if resolved_id is None:
                raise ValueError("SeqRun samples require a resolvable sample_id.")
            if resolved_id in sample_map:
                raise ValueError(f"Duplicate sample_id in SeqRun: {resolved_id}")
            if sample.sample_id != resolved_id:
                raise ValueError(f"SeqRun key/sample_id mismatch: key={resolved_id}, object={sample.sample_id}")
            sample_map[resolved_id] = sample
        return sample_map

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
        return pd.DataFrame([s.sample_meta for s in self.sample_list]) if self.samples else pd.DataFrame()

    @property
    def df(self) -> pd.DataFrame:
        if self.data_df is not None:
            return self.data_df.copy()
        if not self.samples:
            return pd.DataFrame()
        df = pd.concat([s.df for s in self.sample_list], ignore_index=True)
        for col, val in {"run_id": self.run_id, "run_number": self.run_number,
                         "run_name": self.run_name, "run_date": self.run_date}.items():
            if col not in df.columns or df[col].isna().all():
                df[col] = val
        return df

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
        

#%% Versions
"""
v1.0 
- Initial version, accomodating new model/container structure
- Supercedes deprecated parse_countdata module and associated code in bcParse.py
2026-05-26 - Refactored into bcparse package structure
2026-05-26 - Added ensure_runinfo for synthetic compile runinfo
"""
