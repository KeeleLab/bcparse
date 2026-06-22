# xlsx.py
"""
Name:      xlsx.py
Author:     CAG
Version:    1.2
Date:       2026/03/26
Refactored: 2026/05/26
"""

#%% Imports

from dataclasses import dataclass
import os
import re
import pandas as pd
from pathlib import Path
import logging
from typing import Dict, List, Optional, Tuple
from bcparse.qc.identity import make_run_id, normalize_date, normalize_long_df
from bcparse.containers.runinfo import RunInfo

logger = logging.getLogger(__name__)

#%% Classes

class WorkbookParseError(ValueError):
    """
    Raised when an Analysis workbook cannot be parsed into the expected format.
    """
    def __init__(self, message: str, *, filepath: str | Path | None = None) -> None:
        super().__init__(message)
        self.message = str(message)
        self.filepath = None if filepath is None else Path(filepath)

    def with_filepath(self, filepath: str | Path) -> "WorkbookParseError":
        """
        Return a copy of this error bound to a specific workbook path.
        """
        return WorkbookParseError(self.message, filepath=filepath)

    def __str__(self) -> str:
        if self.filepath is None:
            return self.message
        return f"file: {self.filepath}\nreason: {self.message}"

@dataclass
class ParsedAnalysis:
    """
    Container for one parsed Analysis workbook.
    """
    run_id: str
    data_df: pd.DataFrame
    runinfo: RunInfo
    filepath: Path


class XlsxPathManager:
    """
    Initialize the file manager.

    Args:
        wdir: Working directory to search (defaults to the current directory).
        txtfile: Optional text file listing file paths (one per line).
        pattern: Regex pattern to match filenames (defaults to '.*Analysis\\.xlsx$').
    """
    
    def __init__(
        self,
        wdir: Optional[str] = None,
        txtfile: Optional[str] = None,
        pattern: str = r'.*Analysis\.xlsx$'
        ) -> None:

        self.wdir: str = wdir or os.getcwd()
        self.txtfile: Optional[str] = txtfile
        self.pattern: str = pattern
        self.files: List[Path] = self._discover()

    def _discover(self) -> List[Path]:
        """
        Discover input files either from a provided txt file or by scanning the directory.
        """
        files = self._fromtxt(self.txtfile) if self.txtfile else self._fromscan(self.wdir)

        if not files:
            raise FileNotFoundError(
                f"No files matching pattern '{self.pattern}' found in "
                f"{'text file ' + str(self.txtfile) if self.txtfile else 'directory ' + str(self.wdir)}"
            )

        files = self._dedup(files)
        self.files = files
        return files

    def _fromtxt(self, txtfile: str) -> List[Path]:
        """
        Read file paths from a text file and validate their existence.
        """
        with open(txtfile, "r") as f:
            files = [Path(line.strip()) for line in f if line.strip()]

        missing = [p for p in files if not p.is_file()]
        if missing:
            mlist = "\n".join(map(str, missing))
            raise FileNotFoundError(
                f"The following input files were not found:\n{mlist}"
                )

        return files

    def _fromscan(self, directory: str) -> List[Path]:
        """
        Recursively scan the given directory for files matching the pattern.
        """
        directory_path = Path(directory)
        return [
            path
            for path in directory_path.rglob("*")
            if path.is_file()
            and not path.name.startswith("~$")
            and re.match(self.pattern, path.name)
            ]

    def _dedup(self, files: List[Path]) -> List[Path]:
        """
        Deduplicate by exact filename (case-sensitive), keeping the newest file.
        """
        seen: dict[str, Path] = {}
        duplicates: List[tuple[Path, Path]] = []

        for path in files:
            fname = path.name
            if fname in seen:
                newer = max([seen[fname], path], key=lambda p: p.stat().st_ctime)
                older = seen[fname] if newer == path else path
                duplicates.append((older, newer))
                seen[fname] = newer
            else:
                seen[fname] = path

        if duplicates:
            logger.warning("Duplicate filenames found; keeping the newest file:")
            for old, new in duplicates:
                logger.warning(f"  Older: {old}, Newer: {new}")

        return list(seen.values())


class AnalysisWorkbookParser:
    """
    AnalysisWorkbookParser class:
    Parses an Excel analysis file, extracts metadata and sample tables,
    and normalizes column names.
    """

    def __init__(self, filepath: str | Path) -> None:

        self.filepath: Path = Path(filepath)
        try:
            self._load_xlsx()
            self._parse_xlsx()
        except WorkbookParseError as e:
            raise e.with_filepath(self.filepath) from e
        except Exception as e:
            raise WorkbookParseError(
                f"Unexpected {type(e).__name__}: {e}",
                filepath=self.filepath,
            ) from e

    def _load_xlsx(self, **kwargs) -> None:
        """
        Load all sheets from the Excel file into self.sheets.
        """
        #logger.debug(f"Loading Excel file: {self.filepath}")
        self.sheets: Dict[str, pd.DataFrame] = pd.read_excel(
            self.filepath,
            sheet_name=None,  # load all sheets into a dict
            header=None,
            engine="openpyxl",
            **kwargs
            )

    def _parse_xlsx(self) -> None:
        """
        Process the loaded Excel workbook into a standardized long-format DataFrame.
        """
        all_samples: List[pd.DataFrame] = []

        # Obtain runinfo from first sheet
        first_sheet = next(iter(self.sheets.values()))
        self.runinfo = RunInfo.from_table(
            raw_df=self._sheet_to_runinfo_df(first_sheet),
            filepath=str(self.filepath),
        )
        #logger.debug(f"Runinfo: {self.runinfo.run_id}")

        # Process only 'samples' sheets (skip first and non-sample sheets)
        for i, (sname, df) in enumerate(self.sheets.items()):
            if i == 0:
                continue
            if 'samples' not in sname.lower():
                continue

            #logger.debug(f"Parsing: {sname}")
            samples = self._parse_samples(df)
            all_samples.extend(samples)

        self.run_id = make_run_id(
            self.runinfo.run_number,
            self.filepath.name,
            self.filepath,
        )

        self.data_df = self._combine_sample_tables(all_samples)
        self.analysis = ParsedAnalysis(
            run_id=self.run_id,
            data_df=self.data_df.copy(),
            runinfo=self.runinfo,
            filepath=self.filepath,
        )

    def _sheet_to_runinfo_df(self, first_df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert the first sheet of an analysis workbook back into a headered runinfo table.
        """
        if first_df.empty:
            return pd.DataFrame()

        header = first_df.iloc[0].fillna("").astype(str).str.strip()
        runinfo_df = first_df.iloc[1:].reset_index(drop=True).copy()
        runinfo_df.columns = header
        runinfo_df = runinfo_df.loc[:, runinfo_df.columns != ""]
        return runinfo_df

    def _parse_samples(self, df: pd.DataFrame) -> List[pd.DataFrame]:
        """
        Extract multiple side-by-side tables from a 'samples' sheet DataFrame.
        """
        r, c = self._detect_breaks(df)
        
        # Recognizes specific sheet formats:
        valid_breaks = {
            (3, 6), # r=3, c=6: original missing 'short bcs'
            (3, 7), # r=3, c=7: original 7 data columns
            (7, 7), # r=7, c=7: current non-legacy analysis output
            (8, 7), # r=8, c=7: pre-refactor non-legacy analysis output
        }
        if (r, c) not in valid_breaks:
            raise WorkbookParseError(
                "Nonstandard format in samples sheet "
                f"(r={r}, c={c})."
            )

        header_row = df.iloc[r + 1]

        # Define data column blocks
        mask = header_row.notna()
        block_ids = (mask != mask.shift(fill_value=False)).cumsum() * mask
        col_blocks = [
            group.index.tolist()
            for _, group in header_row[mask].groupby(block_ids)
            ]

        tables: List[pd.DataFrame] = []

        for block in col_blocks:
            meta_df = df.iloc[0:r, block]
            meta_dict = self._extract_metadata(meta_df, r)

            block_df = df.iloc[r + 2:, block].copy()
            block_df.columns = df.iloc[r + 1, block]

            # Define above_cutoff and drop spacer rows
            is_blank = block_df.isna().all(axis=1)
            block_df["above_cutoff"] = (is_blank.cumsum() == 0) # this pulls first block, non-derivative and > 1/input
            block_df = block_df[~is_blank].reset_index(drop=True)

            sample_meta = self._resolve_runinfo_sample(meta_dict)
            for k, v in sample_meta.items():
                block_df[k] = v

            block_df["run_number"] = self.runinfo.run_number

            # Add shorts col if missing
            if c == 6:
                block_df["Short barcode?"] = "not checked!"

            # Standardize column names per sample block so each table remains one SeqSamp.
            block_df = block_df.rename(columns={
                "Barcode": "bc_name",
                "Sequence": "bc_seq",
                "Counts": "bc_count",
                "Proportion": "proportion",
                "Hamming Dist to Another Sequence": "ldist_samp_lvl",
                "putative_parent": "ldist_samp_lvl",
                "Index hopping?": "x_idx",
                "Multi-group share": "x_idx",
                "Multi-group?": "x_idx",
                "x_group": "x_idx",
                "Short barcode?": "short_bc"
                })

            if "x_idx" in block_df.columns:
                block_df["x_idx"] = pd.to_numeric(block_df["x_idx"], errors="coerce")

            tables.append(block_df.reset_index(drop=True))

        return tables

    def _combine_sample_tables(self, sample_tables: List[pd.DataFrame]) -> pd.DataFrame:
        """
        Concatenate parsed sample blocks into one normalized long-format dataframe.
        """
        if sample_tables:
            df = pd.concat([table.copy() for table in sample_tables], ignore_index=True)
        else:
            df = pd.DataFrame()

        df["filename"] = self.filepath.name
        df["run_id"] = self.run_id

        df = normalize_long_df(
            df,
            run_number=self.runinfo.run_number,
            run_name=self.runinfo.run_name,
            run_date=self.runinfo.run_date,
            source_path=self.filepath,
        )

        df["run_id"] = df["run_id"].where(df["run_id"].notna(), self.run_id)

        return df

    def _detect_breaks(self, df: pd.DataFrame) -> Tuple[Optional[int], Optional[int]]:
        """
        Detect header break row (r) and determine number of data columns (c).

        Returns:
            Tuple of (r, c) where:
            - r is the row index of the break
            - c is the count of header columns
        """
        mask = df.isna().all(axis=1).to_numpy()
        r = mask.argmax() if mask.any() else None

        if r is not None:
            header_row = df.iloc[r + 1]
            if header_row.isna().any():
                c = header_row.isna().to_numpy().argmax()
            else:
                c = header_row.notna().sum()
        else:
            c = None

        return r, c

    def _extract_metadata(self, meta_block: pd.DataFrame, r: int) -> Dict[str, Optional[str]]:
        """
        Extract sample metadata from a block given the detected break row.
        """
        if r == 3:
            return {
                "samp_group": meta_block.iat[0, 0],
                "samp_name": meta_block.iat[0, 1],
                "samp_date": normalize_date(meta_block.iat[0, 2]),
                "input": meta_block.iat[2, 0],
                }
        elif r == 7:
            return {
                "samp_group": meta_block.iat[3, 1],
                "samp_name": meta_block.iat[4, 1],
                "samp_date": normalize_date(meta_block.iat[5, 1]),
                "input": meta_block.iat[6, 1],
                }
        elif r == 8:
            return {
                "samp_group": meta_block.iat[4, 1],
                "samp_name": meta_block.iat[5, 1],
                "samp_date": normalize_date(meta_block.iat[6, 1]),
                "input": meta_block.iat[7, 1],
                }
        else:
            raise ValueError(f"Unrecognized header break row r={r}")

    def _resolve_runinfo_sample(
        self,
        meta_dict: Dict[str, Optional[str]],
    ) -> Dict[str, object]:
        """
        Recover runinfo-derived sample metadata for a sample block.
        """
        match_df = self.runinfo.sample_rows
        for key in ["samp_group", "samp_name", "samp_date"]:
            match_df = match_df.loc[match_df[key] == meta_dict.get(key)]

        if len(match_df) != 1:
            raise WorkbookParseError(
                "Could not uniquely match sample sheet metadata to runinfo."
            )

        row = match_df.iloc[0]
        return row.to_dict()

#%%
"""
V1.2 
- add normalize_date() as global function from readers, rather than class function
- modified to accomodate new downstream model class container runseries.py
2026-05-26 - Refactored into bcparse package structure
"""
