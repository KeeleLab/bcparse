# io.py

"""
Name:       io.py
Author:     CAG
Version:    1.9
Date:       2026/02/24
"""

# %% Imports

import os
import csv
import gzip
import pandas as pd
import numpy as np
from tqdm import tqdm
from typing import Generator, Tuple
from utils.parse_countdata import parse_countdata
from utils.parse_fastq import parse_fastq_rec
from utils.xlsx_compiler import set_compiler
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

# %% Functions


def read_csv_to_dict(file_path: str):
    """
    Parse two-col csv into dict {'col1':'col1'}
    (Used for primer ref data, k=seq v=name)
    """
    result = {}
    with open(file_path, "r") as csvfile:
        reader = csv.reader(csvfile)
        # Skip the header row
        next(reader)
        for row in reader:
            col1, col2 = row
            result[col1] = col2
    return result


def read_fasta_to_dict(file_path: str):
    """
    Parse a standard fasta into a dict
    (Used for barcode ref data, k=name v=seq)
    """
    sequences = {}
    current_seq = []
    current_header = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                    current_seq = []
                current_header = line[1:]
            else:
                current_seq.append(line)

        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    return sequences


def format_refdata(runinfo_path: str, primer_path: str, barcode_path: str):
    """
    Generate reference runinfo, primer indexes, and barcodes from sources.
    """
    print("Formatting reference data... \n")

    # Load runinfo
    runinfo = pd.read_excel(runinfo_path)

    # Strip whitespace from column names
    runinfo.columns = runinfo.columns.str.strip()

    # Strip whitespace from string data
    runinfo = runinfo.apply(lambda col: col.apply(lambda x: x.strip() if isinstance(x, str) else x))

    # Convert 'Date' column to datetime.date
    runinfo['Date'] = pd.to_datetime(runinfo['Date']).dt.date

    # Replace string 'inf' (any case) with np.inf first
    col = 'Input TOTAL PER BARCODE'
    runinfo[col] = runinfo[col].replace(to_replace=r'(?i)^inf$', value=np.inf, regex=True)

    # Mask for invalid entries
    invalid = (
        ~runinfo[col]
        .astype(str)
        .str.strip()
        .str.match(
            r'^(inf|nan|\d+(\.\d+)?|\.\d+)$', 
            na=False
            ))

    if invalid.any():
        invalid_strs = [str(x) for x in runinfo.loc[invalid, col].unique()]
        raise ValueError(f"Invalid input: {invalid_strs} - use numeric or inf.")

    # Convert column to numeric, errors='raise' should now succeed
    runinfo[col] = pd.to_numeric(runinfo[col], errors='raise')

    # Load primer and barcode reference dicts (assuming functions exist)
    p5_refdict = read_csv_to_dict(primer_path)
    bc_refdict = read_fasta_to_dict(barcode_path)

    return runinfo, p5_refdict, bc_refdict


def stream_fastq(file_path: str) -> Generator[Tuple[str, str, str, str], None, None]:
    """
    Iteratively stream fq or fq.gz entries into memory as 4-line string.
    Shows an indeterminate tqdm progress bar (no pre-counting).
    """
    gz = file_path.endswith(".gz")

    if gz:
        # Buffer the underlying fileobj; gzip itself then reads bigger chunks
        raw = open(file_path, "rb", buffering=1024 * 1024)
        fh = gzip.open(raw, "rt", newline="")
    else:
        fh = open(file_path, "rt", buffering=1024 * 1024, newline="")

    with fh as fastq_file, tqdm(desc="Counting index/barcode pairs", unit="rec", ascii="-=") as pbar:
        it = iter(fastq_file)
        update_every = 1024
        n = 0

        while True:
            h = next(it, None)
            if h is None:
                break
            s = next(it); p = next(it); q = next(it)

            # cheaper than .strip(); we only want to drop newline
            yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")

            n += 1
            if n % update_every == 0:
                pbar.update(update_every)

        # flush remainder
        rem = n % update_every
        if rem:
            pbar.update(rem)


def process_fastq_stream( 
    settings: dict, 
    ):
    """
    Use stream_fastq() and parse_fastq_rec class to collect {idx:{bc:count}}
    """
    # Initialize dict to hold results
    res_pfq = {}

    # Initialize search patterns
    parse_fastq_rec.init_patterns(settings)

    # Let it rip
    for header, seq, plus, qual in stream_fastq(settings["sample_path"]):

        # Parse record to return results.
        res = parse_fastq_rec(
            seq=seq,  # sequence
            qual=qual,  # qstring
            settings=settings,  # run settings
            )

        # Ignore read if barcode is missing or low mean q
        if not hasattr(res, "bc_seq") or res.bc_qual < settings["mean_qual"]:
            continue

        # Process p5 - skip if missing or low mean q
        if not (hasattr(res, "p5_seq") and hasattr(res, "p5_qual")):
            continue

        p5_seq = res.p5_seq if res.p5_qual >= settings["mean_qual"] else None
        if p5_seq is None:
            continue

        # Optionally process p7, fill instead of skip:
        # Ns where q<mask_qual, or all Ns if mean q<min_qual, X if missing
        p7_seq = ""

        if settings["dualindex"]:
            if hasattr(res, "p7_seq") and hasattr(res, "p7_qual"):
                p7_seq = res.p7_seq if res.p7_qual >= settings["mean_qual"] else "N" * settings["tlen_p7"]
            else:
                p7_seq = "X" * settings["tlen_p7"]

        # Set idx_seq as p5::p7, p7 is empty unless dualindex
        idx_seq = p5_seq + p7_seq

        # Add idx subdict to count dict if not present.
        if idx_seq not in res_pfq:
            res_pfq[idx_seq] = {}

        # Add bc to idx subdict if not present.
        if res.bc_seq not in res_pfq[idx_seq]:
            res_pfq[idx_seq][res.bc_seq] = 0

        # Increment counter for {idx:{bc:count}}.
        res_pfq[idx_seq][res.bc_seq] += 1

    return res_pfq


def write_pcd_output(
    runinfo_df: pd.DataFrame,
    settings: dict,
    res_pcd: parse_countdata,
):
    """
    Write the long data as <runname>_concat.csv, and a legacy excel file as <runname>_Analysis.csv
    Formatting is done via parse_countdata instance methods, see countdata_parser.py
    """
    # Pull these from settings
    out_path = settings["out_path"]
    legacy_format = settings["legacy_format"]
    filt_mat_ac = settings["filt_mat_ac"]
    
    def modify_runinfo(runinfo: pd.DataFrame, res_pcd: parse_countdata, missing: list) -> pd.DataFrame:
        """
        Modify runinfo DataFrame with additional analysis metadata
        
        Parameters:
        -----------
        runinfo : pd.DataFrame
            Original runinfo DataFrame
        res_pcd : parse_countdata
            Parsed count data object
        
        Returns:
        --------
        pd.DataFrame
            Modified runinfo DataFrame with additional columns
        """
        # Create a copy to avoid modifying the original DataFrame
        modified_runinfo = runinfo.copy()
        
        # =============== RUNINFO MODIFICATIONS ===============
        # Total Sequences Extracted per Primer (n bc detected w/ valid index)
        modified_runinfo['n_reads'] = modified_runinfo['full_idx'].map(
            {index: sum(barcode.values()) for index, barcode in res_pcd.idx_known.items()}
            )
        
        # Total distinct barcodes
        modified_runinfo['n_bc'] = modified_runinfo['full_idx'].map(
            {index: len(barcodes) for index, barcodes in res_pcd.idx_known.items()}
            )
        
        # Total named barcodes
        modified_runinfo['n_bc_named'] = modified_runinfo['full_idx'].map(
            res_pcd.df[~res_pcd.df["bc_name"].str.contains("Unique", na=False)]
            .groupby("idx_name").size().to_dict()
            )
        
        # Total counts above cutoff
        modified_runinfo['n_reads_ac'] = modified_runinfo['full_idx'].map(
            res_pcd.df[(res_pcd.df["above_cutoff"] == True)]
            .groupby("idx_name")["bc_count"].sum().to_dict()
            )
        
        # Count to input ratio
        modified_runinfo['count-to-input_ratio'] = (
            modified_runinfo['n_reads'] / modified_runinfo['Input TOTAL PER BARCODE']
            )

        # Add column if there are any notes
        if any(res_pcd.idx_notes.values()):
            modified_runinfo['notes'] = modified_runinfo['full_idx'].map(
                {index: '; '.join(notes) for index, notes in res_pcd.idx_notes.items() if notes}
                )

        return modified_runinfo

    # =============== END RUNINFO MODIFICATIONS ===============

    # Create the output directory if it doesn't exist
    os.makedirs(out_path, exist_ok=True)

    # Write long df
    res_pcd.df.to_csv(os.path.join(out_path, f"{res_pcd.rnam}_concat.csv"), index=False)

    # Obtain matrix and table outputs
    group_mats_dict = res_pcd.samp_group_matrices(filt_ac=filt_mat_ac, collapse_to_parent=True)
    group_tabs_dict, _ = res_pcd.samp_group_sidelong_tables(legacy_format=legacy_format,collapse_to_parent=True)
    group_tabs_raw_dict, _ = res_pcd.samp_group_sidelong_tables(legacy_format=legacy_format,collapse_to_parent=False)

    # Get samp_groups from runinfo & output data - find/report missing (no barcodes!)
    expected_groups = runinfo_df["Animal"].dropna().unique()
    groups = res_pcd.df["samp_group"].unique()

    missing = set(expected_groups) - set(groups)

    if missing:
        print(f""" 
        No barcodes found for {missing}!
        """)

    # Check sample groups from runinfo against keys in res dicts
    if not set(group_mats_dict.keys()) == set(groups) or not set(group_tabs_dict.keys()) == set(groups):
        raise ValueError(f"""
        Check concat output for missing data. 
        Details:
        - group_mats_dict keys: {set(group_mats_dict.keys())}
        - group_tabs_dict keys: {set(group_tabs_dict.keys())}
        - expected groups: {set(groups)}

        Previously seen issues here:
        - ran 239M/M2 instead of 239M/M2_dual-index.
        """)

    # Modify runinfo with nested function
    runinfo_df = modify_runinfo(runinfo_df, res_pcd, missing)

    # Open excel writer for xlsx output
    with pd.ExcelWriter(os.path.join(out_path, f"{res_pcd.rnam}_Analysis.xlsx")) as writer:
        
        # Write the runinfo as sheet 1
        runinfo_df.to_excel(writer, sheet_name="runinfo", index=False)

        # For each samp_group...
        for samp_group in groups:
            print(f"Writing {samp_group}...")

            # Write tables - raw
            group_tabs_raw_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} raw", index=False, header=False
                )

            # Write tables - collapsed
            group_tabs_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} samples", index=False, header=False
                )

            # Write matrix
            group_mats_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} matrix", index=True, header=True
                )

        # Add run metadata and arguments to settings_dict and create df
        settings["analysis date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        settings_df = pd.DataFrame(list(settings.items()))

        # Add settings sheet to excel
        settings_df.to_excel(
            writer, sheet_name="Analysis settings", index=False, header=False
            )


def read_base_csv(
    base_csv_path: str,
    ) -> pd.DataFrame:
    """
    Read in and format a previously compiled all.csv file
    """

    return (
        pd.read_csv(base_csv_path)
          .drop(columns=["ldist_group_lvl", "ldist_all"], errors="ignore")
        )


def write_sc_output(
    sc: set_compiler,
    settings: dict,
    ) -> None:
    """
    Format and save final above-cutoff results to an Excel workbook.

        Includes:
        - all_samples sheet (above-cutoff long-form data)
        - one "<samp_group> samples" sheet per group
        - one "<samp_group> matrix" sheet per group (pivoted proportions)
    """
    # Pull these from settings dict
    out_path = settings["out_path"]
    out_prefix = settings["out_prefix"]

    # Build paths
    csv_path = os.path.join(out_path, f"{out_prefix}_all.csv")
    xlsx_path = os.path.join(out_path, f"{out_prefix}_above_cutoff.xlsx")

    # Write full csv
    sc.df.to_csv(csv_path, index=False)

    # Generate per-group above-cutoff matrices
    ac_df = sc.df_ac.copy()
    df_dict, mat_dict = sc.group_mat_dict()

    # Write to Excel
    with pd.ExcelWriter(xlsx_path) as writer:
        ac_df.to_excel(writer, sheet_name="all_samples", index=False)

        if sc.contam_check:
            full_mat, contam_table = sc.contam_report()
            full_mat.to_excel(writer, sheet_name="full_matrix", index=True, header=True)
            contam_table.to_excel(writer, sheet_name="contam_report", index=True, header=True)

        for samp_group in mat_dict:
            print(f"Writing {samp_group} ...")
            df_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} samples", index=False, header=True
                )
            mat_dict[samp_group].to_excel(
                writer, sheet_name=f"{samp_group} matrix", index=True, header=True
                )

#%% Versions
"""
v1.6 
- Added write_sc_output()
- Added validation step for runinfo inputs! inf & numeric okay, others raise error

v1.7
- Modifed parmeter passing from parse and compile methods, affecting multiple functions!

v1.8
- Read in base csv from path

v1.9
- Modified process_fastq_steam() and stream_fastq() for modest speed boost
"""
