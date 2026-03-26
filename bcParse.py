#! /usr/bin/env python3
# bcParse.py

# %% CLI args & logger

import logging
import argparse
from argparse import RawTextHelpFormatter
from dataclasses import asdict, dataclass
from config import stocks, get_settings

ver = "4.0 - 2026.03.26"

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
python bcParse.py --compile --xlsx_path path/to/parent_dir --out_path path/to/output_dir --out_prefix my_project

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
    "--collapse_to_parent",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="Collapse grouped parse outputs to putative parent barcodes (default: True)"
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


# %% Classes

@dataclass
class ParseSettings:
    sample_path: str = ""
    runinfo_path: str = ""
    out_path: str = ""
    stock: str = ""
    dualindex: bool = False
    mean_qual: int = 30
    mismatches: int = 1
    mask: bool = False
    collapse: bool = False
    mask_qual: int = 20
    dist_threshold: int = 1
    filt_mat_ac: bool = False
    legacy_format: bool = False
    collapse_to_parent: bool = True

    @classmethod
    def from_mapping(cls, values: dict) -> "ParseSettings":
        defaults = cls()
        return cls(**{k: values.get(k, getattr(defaults, k)) for k in cls.__annotations__})

    def validate(self) -> None:
        required = ["sample_path", "runinfo_path", "out_path"]
        missing = [field for field in required if not getattr(self, field)]
        if missing:
            raise ValueError(
                "Error: missing args. Call with -h to see full args, or with --gui for graphic interface."
            )

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class CompileSettings:
    xlsx_path: str = ""
    out_path: str = ""
    base_csv_path: str | None = None
    out_prefix: str = ""
    ldist_samp_lvl: bool = False
    ldist_group_lvl: bool = False
    contam_check: bool = False
    dist_threshold: int = 1
    collapse_to_parent: bool = True

    @classmethod
    def from_mapping(cls, values: dict) -> "CompileSettings":
        defaults = cls()
        return cls(**{k: values.get(k, getattr(defaults, k)) for k in cls.__annotations__})

    def validate(self) -> None:
        required = ["xlsx_path", "out_path", "out_prefix"]
        missing = [field for field in required if not getattr(self, field)]
        if missing:
            raise ValueError(
                "Error: missing req'd args. Call with -h to see full args, or with --gui for graphic interface."
            )

    def to_dict(self) -> dict:
        return asdict(self)

# %% Modes

# --- PARSE MODE ---
def run_parse_mode(args):

    # Mode-specific imports
    from utils.parse_fastq import process_fastq_stream
    from utils.seqrun import SeqRun
    
    """
    1. Initialize arguments and run settings. See config.py, utils/io.py, utils/gui.py

    - Required args are sample_path, runinfo_path, out_path, and stock.
    - See config.py for stock options, e.g. "239M"
    - Remaining args have defaults set - see above. Modify if you know what you're doing! 
    """

    # If running with --gui...
    if args.gui:
        from utils.gui import ArgSelectionGUI

        # set variables returned from tkinter object,
        params = ArgSelectionGUI(stocks).run()
        print("GUI params:", params)  

    # otherwise use CLI args.
    else:
        params = ParseSettings(
            sample_path = getattr(args, "sample_path"),
            runinfo_path = getattr(args, "runinfo_path"),
            out_path = getattr(args, "out_path"),
            stock = getattr(args, "stock"),
            dualindex = getattr(args, "dualindex"),
            mean_qual = getattr(args, "mean_qual"),
            mismatches = getattr(args, "mismatches"),
            mask = getattr(args, "mask"),
            collapse = getattr(args, "collapse"),
            mask_qual = getattr(args, "mask_qual"),
            dist_threshold = getattr(args, "dist_threshold"),
            filt_mat_ac = getattr(args, "filt_mat_ac"),
            legacy_format = getattr(args, "legacy_format"),
            collapse_to_parent = getattr(args, "collapse_to_parent"),
        ).to_dict()

    try:
        parse_settings = ParseSettings.from_mapping(params)
        parse_settings.validate()
    except ValueError as e:
        print(str(e))
        return

    params = parse_settings.to_dict()

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
    2. Format reference data and validate runinfo before counting
    """
    runinfo, p5_refdict, bc_refdict, p7_refdict = SeqRun.load_parse_reference_data(
        runinfo_path = settings["runinfo_path"],
        primer_path = settings["primer_path"],
        barcode_path = settings["barcode_path"],
        p7_path = settings["p7_path"] if settings["dualindex"] else None,
    )

    """
    3. Stream fastq to count barcode reads per index

    - Return {idx:{bc:count}}
    - See utils/io.py and utils/fastq_parser.py
    - mask returns N for base q<mask_qual
    - reads drop when mean phred<mean_qual
    - mismatches = # allowed errors
    """

    # Check for sample name collision before counting
    collisions = runinfo.find_idx_sample_id_collisions()
    if collisions:
        lines = [
            "Error: duplicate sample identities found in runinfo.",
            "Distinct idx_name values resolve to the same sample_id:",
            "",
        ]
        for sample_id, idx_names in collisions.items():
            lines.append(f"- {sample_id}")
            lines.append(f"  idx_name(s): {', '.join(idx_names)}")

        raise SystemExit("\n".join(lines))
    else:

        print("\n")

        res_pfq = process_fastq_stream(
            settings=settings, 
            )

        print("\n")

    """
    4. Main analysis: build one SeqRun from {idx_seq:{bc_seq:count}}
    """

    # Build parse-mode helper, then take the core SeqRun from it.
    parse_run = SeqRun.parse_fastq(
        countdict = res_pfq,
        runinfo = runinfo,
        bc_refdict = bc_refdict,
        p5_refdict = p5_refdict,
        p7_refdict = p7_refdict,
        settings = settings,
        )
    seq_run = parse_run.seq_run

    print("\n")

    """
    5. Main output.
    """
    parse_run.write_output(settings=settings, runinfo_df=None)

    # In case of -i run, to examine:
    return res_pfq, seq_run, settings


# --- COMPILE MODE ---
def run_compile_mode(args):

    # Mode-specific imports:
    from utils.parse_xlsx import(
        XlsxPathManager,
        AnalysisWorkbookParser,
        )

    from utils.runseries import RunSeries
    
    # Define expected parameters
    # If running with --gui...
    if args.gui:
        from utils.gui import CompileArgSelectionGUI

        # set variables returned from tkinter object,
        params = CompileArgSelectionGUI().run()
        print("GUI params:", params)

    # otherwise use CLI args.
    else:
        params = CompileSettings(
            xlsx_path = getattr(args, "xlsx_path"),
            out_path = getattr(args, "out_path"),
            base_csv_path = getattr(args, "base_csv_path"),
            out_prefix = getattr(args, "out_prefix"),
            ldist_samp_lvl = getattr(args, "ldist_samp_lvl"),
            ldist_group_lvl = getattr(args, "ldist_group_lvl"),
            contam_check = getattr(args, "contam_check"),
            dist_threshold = getattr(args, "dist_threshold"),
            collapse_to_parent = getattr(args, "collapse_to_parent"),
        ).to_dict()

    try:
        compile_settings = CompileSettings.from_mapping(params)
        compile_settings.validate()
    except ValueError as e:
        print(str(e))
        return

    params = compile_settings.to_dict()

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
    xpm = XlsxPathManager(wdir=settings["xlsx_path"])
    print(f"Discovered {len(xpm.files)} Analysis workbook(s).")

    # Optional base csv in
    base_df = RunSeries.read_base_csv(settings["base_csv_path"]) if settings["base_csv_path"] else None
    if base_df is not None:
        print(f"Loaded base csv with {len(base_df)} row(s).")

    # Iterate paths and obtain list of parsed xlsx data
    res_px = []

    for i, f in enumerate(xpm.files, start=1):
        print(f"Parsing workbook {i}/{len(xpm.files)}: {f.name}")
        try:
            px = AnalysisWorkbookParser(filepath=f)
            if not px.data_df.empty:
                res_px.append(px.analysis)
        except Exception as e:
            logger.error(f"Failed processing {f}: {e}")

    print(f"Parsed {len(res_px)} non-empty analysis workbook(s).")

    # Compile set of parsed data to final output
    print("Building compiled series...")
    run_series = RunSeries.build_from_parsed_files(
        parsed_files = res_px,
        base_df = base_df,
        settings = settings,
        )

    print("Optional QC and output writing...")
    run_series.write_output(settings=settings)
    print("Compile complete.")

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
        level = logging.DEBUG,  # or INFO
        format = "%(asctime)s %(levelname)s %(name)s %(message)s",
        handlers = [logging.StreamHandler()],
        force = True,
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
        res_pfq, seq_run, settings = _ret

        # Optional: print a short confirmation
        print(
            "\nRun with -i (interactive)? Python env will retain:"
            "\n  res_pfq   → raw FASTQ parse results (dict)"
            "\n  seq_run   → parse-mode SeqRun"
            "\n  settings  → run settings dictionary\n"
        )

# %% Changes & versions
"""
****************
* Nice to have *
****************

- Write arg for xlsx path txtfile list
- Streamline running multiple parse runs?

- Brandon wants a QC app... mode?
- Develop this software into a more comprehensive analysis suite?
    - "live" mode, build a list of `bcsamp` objects, standardized figs, etc... 

************
* Versions *
************

v4.0
0. Major overhaul of entire codebase, moving to class-based containers for runinfo, seqsamp, seqrun, and runseries. 
    - All files drastically affected.
1. Dual-index X/INT config
2. Changed parent assignment logic for uniques: named > unique to highest named) (qc_ops.py 2.1->2.2)
3. To qc_ops
    - added informative suffix to parent name
    - big: core_bc search for unresolved uniques in flag_putative_parents
4. Dropped matrix backend outright, no need to carry forward
5. Change x_group to x_idx, and make it numeric for easier thresholding/idx hopping checks
6. Pre-flight parse mode sample name collision check, with informative error message
7. Optional collapse to parent in both parse and compile modes, default True. See settings and qc_ops.py 2.2 for details.
5. Resetting my comment log - previous edits can be found in deprecated branch bcParse_v3.42

"""