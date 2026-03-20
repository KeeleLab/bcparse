#! /usr/bin/env python3
# bcParse.py

# %% CLI args & logger

import logging
import argparse
from argparse import RawTextHelpFormatter
from config import stocks, get_settings

ver = "3.43 - 2026.03.20"

description = f"""
Name:       bcParse.py
Author:     CAG
Version:    {ver}

This program runs two distinct processes:

"Parse": 
Given a short-read fastq file and associated runinfo csv, parses data into a 
formatted excel workbook and long-format csv file. 

"Compile":
Given a directory containing "*Analysis.xlsx" files (produced by Parse mode),
compile all analysis files into a single csv and corresponding excel file

Usage:
Run the script in CLI mode by default or use the --gui flag to launch the Tkinter GUI.
Use the --compile flag (--gui compatible) to run compile mode. 

Example CLI command:
python bcParse.py --sample_path path/to/your.fastq --runinfo_path path/to/runinfo.xlsx --out_path path/to/wherever --stock 239M --dualindex
python bcParse.py --compile --xlsx_path path/to/parent_dir --compile_out_csv path/to/filename.csv --compile_out_xlsx path/to/filename.xlsx

Example GUI command:
python bcParse.py --gui # will open a dialogue window for moused selection
python bcParse.py --gui --compile # to open a compiler-specific dialogue window
"""

# Setup logger
logger = logging.getLogger(__name__)

# Create the parser
parser = argparse.ArgumentParser(
    description=description, formatter_class=RawTextHelpFormatter
    )

# --- PARSE MODE ARGS ---
parser.add_argument(
    "--sample_path", type = str,
    help="path/to/sample.fq"
    )

parser.add_argument(
    "--runinfo_path", type = str,
    help="path/to/runinfo.xlsx"
    )

parser.add_argument(
    "--stock", choices=stocks,
    help="Declare stock-specific run mode"
    )

parser.add_argument(
    "--dualindex", action="store_true",
    help="Run in dual-index mode"
    )

parser.add_argument(
    "--mean_qual", type=int, default=30,
    help="Minimum allowed mean phred per read/seq extract"
    )

parser.add_argument(
    "--mismatches", type=int, default=1,
    help="Number of allowed target mismatches for fastq parsing"
    )

parser.add_argument(
    "--mask", action="store_true",
    help="Mask low quality bases"
    )

parser.add_argument(
    "--mask_qual", type=int, default=20,
    help="Mask low quality bases (default: 20)"
    )

parser.add_argument(
    "--collapse", action="store_true",
    help="Collapse single-N ambiguous reads"
    )

parser.add_argument(
    "--legacy_format", action="store_true",
    help="Return output in legacy format"
    )

parser.add_argument(
    "--filt_mat_ac", action="store_true",
    help="Filter matrix for above-cutoff only"
    )

parser.add_argument(
    "--gui", action="store_true",
    help="Use GUI for file selection"
    )

# --- COMPILE MODE ARGS ---
parser.add_argument(
    "--compile", action="store_true",
    help="Run compile protocol instead of parsing"
    )

parser.add_argument(
    "--xlsx_path", type=str,
    help="path/to/dir with Analysis.xlsx files"
    )

parser.add_argument(
    "--base_csv_path", type=str, default=None,
    help="Optional: provide a path to a previously-compiled *all.csv to coalesce against new compiled analyses"
)

parser.add_argument(
    "--out_prefix", type=str,
    help="Output prefix for compile mode"
    )

parser.add_argument(
    "--ldist_samp_lvl", action="store_true",
    help="Rerun ldist checks at sample-level"
    )

parser.add_argument(
    "--ldist_group_lvl", action="store_true",
    help="Run ldist checks on compiled groups"
    )

parser.add_argument(
    "--contam_check", action="store_true",
    help="Run contamination check & return extra sheet"
    )

# --- COMMON ARGS ---
parser.add_argument(
    "--out_path", type = str,
    help="path/to/output_location"
    )

parser.add_argument(
    "--dist_threshold", type = int, default = 1,
    help="Distance cutoff for parent checks"
    )

# %% Modes

# --- PARSE MODE ---
def run_parse_mode(args):

    # Mode-specific imports
    from utils.io import (
        read_csv_to_dict,
        format_refdata,
        process_fastq_stream,
        write_pcd_output,
        )
    from utils.parse_countdata import parse_countdata
    
    """
    1. Initialize arguments and run settings. See config.py, utils/io.py, utils/gui.py

    - Required args are sample_path, runinfo_path, out_path, and stock.
    - See config.py for stock options, e.g. "239M"
    - Remaining args have defaults set - see above. Modify if you know what you're doing! 
    """

    FIELD_NAMES = [
        "sample_path",
        "runinfo_path",
        "out_path",
        "stock",
        "dualindex",
        "mean_qual",
        "mismatches",
        "mask",
        "collapse",
        "mask_qual",
        "dist_threshold",
        "filt_mat_ac",
        "legacy_format",
        ]
    
    # If running with --gui...
    if args.gui:
        from utils.gui import ArgSelectionGUI

        # set variables returned from tkinter object,
        values = ArgSelectionGUI(stocks).run()
        params = dict(zip(FIELD_NAMES, values))
        print("GUI params:", params)  

    # otherwise use CLI args.
    else:
        params = {k: getattr(args, k) for k in FIELD_NAMES}

    if not all([params["sample_path"], params["runinfo_path"], params["out_path"]]):
        print(
            "Error: missing args. Call with -h to see full args, or with --gui for graphic interface."
            )
        return

    # Show off...
    print(
    f"""

    ▌     ▄▖▄▖▄▖▄▖▄▖
    ▛▌▛▘  ▙▌▌▌▙▘▚ ▙▖
    ▙▌▙▖▗ ▌ ▛▌▌▌▄▌▙▖
    {ver}
    Keele lab, ACVP, FNLCR

    """
    )

    # Get stock-defined run settings - add new stocks and settings in config.py!
    settings = get_settings(
        params["stock"], 
        params["dualindex"]
        )
    
    # Add CLI/GUI args to settings dict:
    settings.update(params)
    
    settings.update({
        "software_version": ver,
    })

    # Write to terminal
    print("\n".join(f"{key}: {value}" for key, value in settings.items()) + "\n")

    """
    2. Stream fastq to count barcodes reads per index
    
    - Return {idx:{bc:count}}
    - See utils/io.py and utils/fastq_parser.py
    - mask returns N for base q<mask_qual
    - reads drop when mean phred<mean_qual
    - mismatches = # allowed errors
    """
    res_pfq = process_fastq_stream(
        settings=settings, 
        )

    print("\n")

    """
    3. Main analysis: parse {idx_seq:{bc_seq:count}} to multiattr obj res_pcd
    
    - Incl. final processed 'df' & class methods for output 
    - See utils/countdata_parser.py 
    """
    # Format reference data for parse_countdata()
    runinfo, p5_refdict, bc_refdict = format_refdata(
        settings["runinfo_path"], 
        settings["primer_path"], 
        settings["barcode_path"]
        )

    # adjust for dual index runs:
    if settings["dualindex"]:
        p7_refdict = read_csv_to_dict(settings["p7_path"])
        runinfo["full_idx"] = runinfo["Barcodes"] + "_" + runinfo["(F Barcode)"]
    else:
        # Default to p5-only settings
        p7_refdict = None
        runinfo["full_idx"] = runinfo["Barcodes"]

    # Run and hold instance as res_pcd:
    res_pcd = parse_countdata(
        countdict=res_pfq,
        runinfo=runinfo,
        bc_refdict=bc_refdict,
        p5_refdict=p5_refdict,
        p7_refdict=p7_refdict,
        settings=settings,
        )

    print("\n")

    """ 
    4. Main output - see io.py & utils/countdata_parser.py. 
    """
    write_pcd_output(
        runinfo_df=runinfo,
        settings=settings,
        res_pcd=res_pcd,
        )

    # In case of -i run, to examine:
    return res_pfq, res_pcd, settings


# --- COMPILE MODE ---
def run_compile_mode(args):

    # Mode-specific imports:
    import gc
    from pathlib import Path
    from typing import List

    from utils.parse_xlsx import(
        xlsx_pathmanager,
        parse_xlsx
        )

    from utils.xlsx_compiler import set_compiler
    from utils.io import (
        read_base_csv,
        write_sc_output
        )
    
    # Define expected parameters
    FIELD_NAMES = [
        "xlsx_path",
        "out_path",
        "base_csv_path",
        "out_prefix",
        "ldist_samp_lvl",
        "ldist_group_lvl",
        "contam_check",
        "dist_threshold",
        ]

    # If running with --gui...
    if args.gui:
        from utils.gui import CompileArgSelectionGUI

        # set variables returned from tkinter object,
        values = CompileArgSelectionGUI().run()
        params = dict(zip(FIELD_NAMES, values))
        print("GUI params:", params)

    # otherwise use CLI args.
    else:
        params = {k: getattr(args, k) for k in FIELD_NAMES}

    if not all([params["xlsx_path"], params["out_path"], params["out_prefix"]]):
        print(
            "Error: missing req'd args. Call with -h to see full args, or with --gui for graphic interface."
        )
        return

    # Collect settings dict
    settings = params
    settings["core_bc"] = None   # For now, until I figure out how to detect stocks and do core search...
    
    # Show off...
    print(
    f"""

    ▌     ▄▖▄▖▖  ▖▄▖▄▖▖ ▄▖
    ▛▌▛▘  ▌ ▌▌▛▖▞▌▙▌▐ ▌ ▙▖
    ▙▌▙▖▗ ▙▖▙▌▌▝ ▌▌ ▟▖▙▖▙▖
    {ver}
    Keele lab, ACVP, FNLCR

    """
    )

    # Set up paths to Analysis files
    xpm = xlsx_pathmanager(wdir=settings["xlsx_path"])

    # Optional base csv in
    base_df = read_base_csv(settings["base_csv_path"]) if settings["base_csv_path"] else None

    # Iterate paths and obtain list of parsed xlsx data
    res_px = []

    for f in xpm.files:
        try:
            px = parse_xlsx(filepath=f)
            if not px.parsed_data.empty:
                res_px.append(px.parsed_data)
        except Exception as e:
            logger.error(f"Failed processing {f}: {e}")

    # Compile set of parsed data to final output
    sc = set_compiler(
        parsed_files = res_px,
        base_df=base_df,
        settings = settings,
        )

    write_sc_output(
        sc = sc,
        settings = settings,
        )

# %% Main process function

def main():
    """
    Entry point for bcParse.py

    Parses CLI args and runs either Parse or Compile mode.
    Returns key objects so they can be accessed interactively
    (e.g., with `python -i bcParse.py ...`).
    """
    # Configure logger
    logging.basicConfig(
        level=logging.DEBUG,  # or INFO
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        handlers=[logging.StreamHandler()],
        force=True,
        )

    # Parse command-line arguments once
    args = parser.parse_args()

    # Run the correct mode and return its result
    if args.compile:
        return run_compile_mode(args)
    else:
        return run_parse_mode(args)


if __name__ == "__main__":
    # Capture the return from main()
    _ret = main()

    # Expose objects for -i interactive sessions
    if _ret is not None:
        # Unpack exactly in the order returned by run_parse_mode()
        res_pfq, res_pcd, settings = _ret

        # Optional: print a short confirmation
        print(
            "\nRun with -i (interactive)? Python env will retain:"
            "\n  res_pfq   → raw FASTQ parse results (dict)"
            "\n  res_pcd   → parsed count data (parse_countdata instance)"
            "\n  settings  → run settings dictionary\n"
        )

# %% Changes & versions
""" 
**********************
* Requested/Critical *
**********************

- Add index-hopping QC check. 
    - x_idx, and less than 0.3% of the proportion found in a higher sample in the run

- Duplicate name repair/add warning in parse mode 
    - For legacy runs after this, will need to enumerate duplicate sample names...

****************
* Nice to have *
****************

- Write arg for xlsx path txtfile list

- speed fastq_parse by counting unique reads first?
    - pro: Do regex/trim ops once per read set
    - con: lose ability to do per-read nt qc
    - Maybe make that optional based on 'mask' param/if not mask then pile reads?

*************
* Rainy day *
*************

- Brandon wants a QC app... mode?

    - single class to represent a complete sample, `bcsamp`
        - parse_countdata would build these from fastq/runinfo
        - parse_xlsx would build these from xlsx
        - If that could be a standard format, then a lot of operations could be recycled as class methods...

    - Develop this software into a more comprehensive analysis suite?
        - "live" mode, build a list of `bcsamp` objects, standardized figs, etc... 

************
* Versions *
************

v3.43
1. Dual-index X/INT config (config.py 1.5->1.6, parse_fastq.py 1.4->1.5)
2. Changed parent assignment logic for uniques: named > unique to highest named) (qc_ops.py 2.1->2.2)
3. To qc_ops
    - added informative suffix to parent name (qc_ops.py 2.1 -> 2.2 -> 2.3)
    - moved functions from parsefastq to qc_ops
    - big: core_bc search for unresolved uniques in flag_putative_parents
4. Dropped matrix backend outright, no need to carry forward

v3.42
1. Add M/M2+Spike configs & ref data

v3.41
1. Normalize column datatypes in set_compiler._combine_with_base() to harden matrix pivot

v3.4
1. For analysis matrix output, modified to exclude above-cutoff children (omit second "block" in sidelong tables)
2. Added option to include a previously-compiled base 'all.csv' in compile mode
    - Behavior is to drop overlapping identical sample/runs from base and default to the new analysis
3. Fixed logger setup to be global (only I care about this)

v3.32
1. Fixed a major bug in combine_known_idx: 
    - dist 1 between expected indexes (e.g. P5.40, P5.60) caused failed assignment! 

v3.31
1. All I did was change primer index names from "VPX.P5/7.#" to "P5/7.#"

v3.3
1. Removed default masking behavior and separated mask/collapse to seperate params
2. precompile regex patterns for fastq parsing - >100% speed bump!
3. BREAKING - integrated bk_tree parent search, make sure to rerun _third_party_packages.py to get rapidfuzz
4. Dropped Levenshtein package in leu of rapidfuzz.levenshtein across all usage
5. Added dist_backend param, can run either matrix or bk parent-finding
6. Made dist_threshold adjustable and made parent dist reporting dynamic (qc_ops/get_filtered_hits, parse_countdata/legacy format)
8. Brought back dualindex as T/F param
9. Changed param handling significantly, pass 'settings' dict to classes.
10. Edited GUI to accomodate new parameters
11. Included -i returns to examine objects post-run if desired
12. Changed input ratio to total counts instead of above-cutoff.
13. Added BFP9 spike to M/M2 references
14. Created a combined M_plus_M2 reference
15. Added new P7 index sequences

v3.2
1. Repaired broken above-cutoff proportion normalization (was missing samp_group grouping factor)
2. Repaired contam_check missing flag

v3.1
1. Combined common methods across parse_countdata and set_compiler to utils.qc_ops
2. Modified compiler output naming variables
3. Added optional qc steps for compiled data incl. contam_check/full_mat
4. Inf input treated as 10 for cutoff determination, non-numeric throws error
5. Added SHIV Ms to config

! v3.0 !
1. Added compiler mode w/ modules parse_xlsx.py and xlsx_compiler.py
2. Changed naming conventions of earlier modules for clarity. 

v2.53
1. In cases where no named barcodes above cutoff are found, instead use nominal 1/input (min theoretical named) for uniques
2. Conditionally add per-index 'notes' column to runinfo
    - "ambig_children collapsed"
    - "no named above_cutoff"
3. Drop non-collapsed ambig uniques from results

v2.52
1. Updated config to run all "M" stocks with dual-index
2. Changed fastq streaming behavior to report count rather than progress bar - saves a lot of time

v2.51
1. Moved unique/above cutoff to rowgroup 0 for excel sample sheets  
2. Recalculate rowgroup 0 props on ouput sample sheets to total 1

v2.50
1. Added modifiable quality params to adjust minimum req'd phred scores
    - mean_qual, mask_qual, and mismatches, passed through argparse or GUI
2. Changed date formats
3. Sort matrix output by date
4. New default run includes fuzzy matching and lowQ base masking
5. Changed argument handling across accessory scripts, all were reversioned
6. Substantially modified GUI
7. (Non-python) ./_misc_/MacroTemplate.xlsm included to recalculate proportions post-run.

v2.43
1. added samp_date to matrix output
2. modified i/o error message in response to a single-vs-dual index issue w/ identical sample names.
3. modified to handle low input samples such that "empty" indexes with no barcodes don't break pipe. 

v2.42
1. Add columns to output runinfo:
    - n_reads, n_bc, n_reads_ac, n_bc_named, count-to-input_ratio (see io.py/write_output())
2. Modify output matrix:
    - Add sequencing input row to header
    - move bc_name and bc_seq to left

v2.41
1. Made output compatible with R-based PCRU concatination app (samp_date format problem)
2. Added 239X/INT, NIRM, OptM, and dGY stocks to config
3. Strip invisible whitespace from runinfo
4. Error check for known indexes - suggest runinfo format issue
5. Stress-test vs. large NIRM run

v2.40
1. Much time spent on making pyright happy (explicit type declarations, etc) - no functional changes

v2.30
1. removed --dualindex as argument and defined as per-stock setting in config.py (incl updates to gui and args)
2. added 239M/M2_dualindex as stock settings to accomodate 1.
3. BIG: for 239M/M2_dualindex, I extract VPX to VPR as barcode, whatever length, and fill to tlen_bc with vpx_end
        - This is to match previous versions 'short' barcodes
        - I think we should be reporting the actual extracted barcode instead, but whatever!
4. fastq_parser.py is fundamentally different now, incl. new match_extract_dual_targ()

v2.20:
1. parse_fastq_rec returned as instance instead of dict repr
2. p5/p7_seq and p5/p7_qual are held as independent attrs during parse_fastq_rec init, joined (or not) in process_fastq_stream()
3. bases < q20 are optionally substituted with N, resulting in idx and bc_seqs showing N, which can be leveraged for finer-grain analyses
4. --mask, --filt_mat_ac, --legacy_format added as optional boolean params, default False
5. GUI modified to indicate new params
"""
