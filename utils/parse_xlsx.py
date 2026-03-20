# parse_xlsx.py

"""
Name:       parse_xlsx.py
Author:     CAG
Version:    1.1
Date:       2026/02/05
"""

#%% Imports

import os
import re
import pandas as pd
from pathlib import Path
import logging
from typing import Dict, List, Optional, Tuple
from utils.qc_ops import normalize_date

logger = logging.getLogger(__name__)

#%% Classes

class xlsx_pathmanager:
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


class parse_xlsx:
    """
    parse_xlsx class:
    Parses an Excel analysis file, extracts metadata and sample tables,
    and normalizes column names.
    """

    def __init__(self, filepath: str | Path) -> None:

        self.filepath: Path = Path(filepath)
        self._load_xlsx()
        self._parse_xlsx()

    def _load_xlsx(self, **kwargs) -> None:
        """
        Load all sheets from the Excel file into self.sheets.
        """
        logger.debug(f"Loading Excel file: {self.filepath}")
        self.sheets: Dict[str, pd.DataFrame] = pd.read_excel(
            self.filepath,
            sheet_name=None,  # load all sheets into a dict
            header=None,
            engine="openpyxl",
            **kwargs
            )

    def _parse_xlsx(self) -> None:
        """
        Process the loaded Excel workbook into a standardized parsed_data DataFrame.
        """
        all_samples: List[pd.DataFrame] = []

        # Obtain runinfo from first sheet
        runinfo = self._parse_runinfo(next(iter(self.sheets.values())))
        logger.debug(f"Runinfo meta: {runinfo}")

        # Process only 'samples' sheets (skip first and non-sample sheets)
        for i, (sname, df) in enumerate(self.sheets.items()):
            if i == 0:
                continue
            if 'samples' not in sname.lower():
                continue

            logger.debug(f"Parsing: {sname}")
            samples = self._parse_samples(df)
            all_samples.extend(samples)

        # Define final DataFrame for this xlsx
        self.parsed_data: pd.DataFrame = pd.concat(all_samples, ignore_index=True)
        self.parsed_data["filepath"] = self.filepath
        self.parsed_data["filename"] = self.filepath.name

        for k, v in runinfo.items():
            self.parsed_data[k] = v

    def _parse_runinfo(self, first_df: pd.DataFrame) -> Dict[str, Optional[str]]:
        """
        Extract run metadata from the first sheet of the workbook.
        """
        col0 = first_df.iloc[:, 0].astype(str).str.strip().str.lower()

        # Map input labels to output keys
        label_map = {
            "run number": "run_number",
            "run name": "run_name",
            "date": "run_date"
            }

        result: Dict[str, Optional[str]] = {}
        for label, out_key in label_map.items():
            match_idx = col0[col0 == label].index
            if not match_idx.empty:
                val = first_df.iat[match_idx[0] + 1, 0]
                if label == "date":
                    val = normalize_date(val)
                result[out_key] = val
            else:
                result[out_key] = None

        return result

    def _parse_samples(self, df: pd.DataFrame) -> List[pd.DataFrame]:
        """
        Extract multiple side-by-side tables from a 'samples' sheet DataFrame.

        Recognizes specific sheet formats:
        - r=3, c=6: original missing 'short bcs'
        - r=3, c=7: original 7 data columns
        - r=8, c=7: new format with metadata in rows 0-7
        """
        r, c = self._detect_breaks(df)
        if r not in (3, 8) or c not in (6, 7):
            logger.warning(f"Break detection failed (r={r}, c={c}); skipping sheet.")
            return []

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

            # Add metadata columns
            for k, v in meta_dict.items():
                block_df[k] = v

            tables.append(block_df)

        if not tables:
            return []

        sheet_df = pd.concat(tables, ignore_index=True)

        # Add shorts col if missing
        if c == 6:
            sheet_df["Short barcode?"] = "not checked!"

        # Standardize column names
        sheet_df = sheet_df.rename(columns={
            "Barcode": "bc_name",
            "Sequence": "bc_seq",
            "Counts": "bc_count",
            "Proportion": "proportion",
            "Hamming Dist to Another Sequence": "ldist_samp_lvl",
            "putative_parent": "ldist_samp_lvl",
            "Index hopping?": "multi_idx",
            "Multi-group?": "multi_idx",
            "x_group": "multi_idx",
            "Short barcode?": "short_bc"
            })

        return [sheet_df]

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
        elif r == 8:
            return {
                "samp_group": meta_block.iat[4, 1],
                "samp_name": meta_block.iat[5, 1],
                "samp_date": normalize_date(meta_block.iat[6, 1]),
                "input": meta_block.iat[7, 1],
                }
        else:
            raise ValueError(f"Unrecognized header break row r={r}")

#%%
"""
V1.1 - add normalize_date() as global function from utils.io, rather than class function
"""
