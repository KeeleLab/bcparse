#! /usr/bin/env python3
# __main__.py - main entry point for bcparse package

# %% CLI args & logger

import argparse
import logging
from argparse import RawTextHelpFormatter
from pathlib import Path

# Allow direct full-path execution from outside the repo:
#   python /path/to/repo/bcparse/__main__.py --gui
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bcparse.config import get_settings, stocks
from bcparse.settings import CompileSettings, ParseSettings

ver = "4.2.0 - 2026.05.27"

description = f"""
Name:       bcparse
Author:     CAG
Version:    {ver}

Barcode parsing and compilation for short-read sequencing data.

Two modes:
  Parse:    fastq + runinfo  ->  formatted Excel workbook + long-format csv
  Compile:  directory of *Analysis.xlsx files, optional base long-format csv  ->  single compiled csv + xlsx

------------------------------------------------------------------------------
Installation
------------------------------------------------------------------------------
Install once, then run `bcparse` from any directory.

  uv tool install /path/to/bcParse_v{ver.split()[0]}       (recommended, uv)
  pipx install    /path/to/bcParse_v{ver.split()[0]}       (recommended, pipx)
  pip install     /path/to/bcParse_v{ver.split()[0]}       (standard pip)

The first two create an isolated environment for bcparse and put the command
on your PATH. Plain pip installs into the current environment.

To update, reinstall with --force (uv, pipx) or --upgrade (pip).

------------------------------------------------------------------------------
Usage
------------------------------------------------------------------------------
After install, `bcparse` works from any working directory. The launch
directory only affects where GUI file dialogs open by default.

GUI:
  bcparse --gui                 # parse mode
  bcparse --gui --compile       # compile mode

Parse (CLI):
  bcparse --sample_path path/to/sample.fastq \\
          --runinfo_path path/to/runinfo.xlsx \\
          --out_path path/to/output \\
          --stock 239M --dualindex

Compile (CLI):
  bcparse --compile \\
          --xlsx_path path/to/analysis_dir \\
          --out_path path/to/output \\
          --out_prefix my_project

See --help for the full argument list and stock options.

------------------------------------------------------------------------------
Running without installing
------------------------------------------------------------------------------
If dependencies are already available in the active Python environment,
bcparse can also be invoked directly from the source tree without installing:

  python -m bcparse --gui                            (from project root)
  python /path/to/bcParse_v{ver.split()[0]}/bcparse/__main__.py --gui   (from anywhere)

Both forms accept the same arguments as the installed `bcparse` command.
Useful for quick testing of a source checkout, but installation is preferred
for normal use.

------------------------------------------------------------------------------
Developer setup
------------------------------------------------------------------------------
For working on bcparse itself, use an editable install so code edits are
picked up without reinstalling:

  uv tool install --editable .      (or)
  pipx install   --editable .       (or)
  pip install    -e .

For a fully locked dev environment (matches uv.lock, includes test/lint deps):

  cd bcParse_v{ver.split()[0]}
  uv sync
  uv run bcparse --gui
  uv run pytest
  uv run ruff check .

"""

# Setup logger
logger = logging.getLogger(__name__)

# Create the parser
parser = argparse.ArgumentParser(
    prog="python -m bcparse",
    description=description,
    formatter_class=RawTextHelpFormatter,
)

# --- PARSE MODE ARGS ---
parser.add_argument("--sample_path", type=str, help="path/to/sample.fq")

parser.add_argument("--runinfo_path", type=str, help="path/to/runinfo.xlsx")

parser.add_argument("--stock", choices=stocks, help="Declare stock-specific run mode")

parser.add_argument("--dualindex", action="store_true", help="Run in dual-index mode")

parser.add_argument(
    "--mean_qual",
    type=int,
    default=30,
    help="Minimum allowed mean phred per read/seq extract",
)

parser.add_argument(
    "--mismatches",
    type=int,
    default=1,
    help="Number of allowed target mismatches for fastq parsing",
)

parser.add_argument("--mask", action="store_true", help="Mask low quality bases")

parser.add_argument(
    "--mask_qual", type=int, default=20, help="Mask low quality bases (default: 20)"
)

parser.add_argument(
    "--collapse_ambig", action="store_true", help="Collapse single-N ambiguous reads"
)

parser.add_argument(
    "--legacy_format", action="store_true", help="Return output in legacy format"
)

parser.add_argument(
    "--filt_mat_ac", action="store_true", help="Filter matrix for above-cutoff only"
)

parser.add_argument(
    "--collapse_to_parent",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="Collapse grouped parse outputs to putative parent barcodes (default: True)",
)

parser.add_argument(
    "--append_spike_ref",
    action="store_true",
    help="Append Spike_reference.fa entries to the selected barcode reference in parse mode",
)

parser.add_argument("--gui", action="store_true", help="Use GUI for file selection")

# --- COMPILE MODE ARGS ---
parser.add_argument(
    "--compile", action="store_true", help="Run compile protocol instead of parsing"
)

parser.add_argument(
    "--xlsx_path", type=str, help="path/to/dir with Analysis.xlsx files"
)

parser.add_argument(
    "--base_csv_path",
    type=str,
    default=None,
    help="Optional: provide a path to a previously-compiled *all.csv to coalesce against new compiled analyses",
)

parser.add_argument("--out_prefix", type=str, help="Output prefix for compile mode")

parser.add_argument(
    "--ldist_samp_lvl", action="store_true", help="Rerun ldist checks at sample-level"
)

parser.add_argument(
    "--ldist_group_lvl", action="store_true", help="Run ldist checks on compiled groups"
)

parser.add_argument(
    "--contam_check",
    action="store_true",
    help="Run contamination check & return extra sheet",
)

# --- COMMON ARGS ---
parser.add_argument("--out_path", type=str, help="path/to/output_location")

parser.add_argument(
    "--dist_threshold", type=int, default=1, help="Distance cutoff for parent checks"
)


# %% Modes


# --- PARSE MODE ---
def run_parse_mode(args):

    # Mode-specific imports
    from bcparse.emit.workbooks import write_parse_workbook
    from bcparse.ingest.countdict import CountDictBuilder, load_parse_reference_data
    from bcparse.ingest.fastq import stream_to_countdict

    """
    1. Initialize arguments and run settings. See bcparse/config.py, bcparse/ingest/fastq.py, bcparse/gui/

    - Required args are sample_path, runinfo_path, out_path, and stock.
    - See bcparse/config.py for stock options, e.g. "239M"
    - Remaining args have defaults set - see above. Modify if you know what you're doing!
    """

    # If running with --gui...
    if args.gui:
        from bcparse.gui.parse import ArgSelectionGUI

        # set variables returned from tkinter object,
        params = ArgSelectionGUI(stocks).run()
        print("GUI params:", params)

    # otherwise use CLI args.
    else:
        params = ParseSettings(
            sample_path=getattr(args, "sample_path"),
            runinfo_path=getattr(args, "runinfo_path"),
            out_path=getattr(args, "out_path"),
            stock=getattr(args, "stock"),
            dualindex=getattr(args, "dualindex"),
            mean_qual=getattr(args, "mean_qual"),
            mismatches=getattr(args, "mismatches"),
            mask=getattr(args, "mask"),
            collapse_ambig=getattr(args, "collapse_ambig"),
            mask_qual=getattr(args, "mask_qual"),
            dist_threshold=getattr(args, "dist_threshold"),
            filt_mat_ac=getattr(args, "filt_mat_ac"),
            legacy_format=getattr(args, "legacy_format"),
            collapse_to_parent=getattr(args, "collapse_to_parent"),
            append_spike_ref=getattr(args, "append_spike_ref"),
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

    # Get stock-defined run settings - add new stocks and settings in bcparse/config.py!
    settings = get_settings(
        params["stock"],
        params["dualindex"],
        append_spike_ref=params["append_spike_ref"],
    )

    # Add CLI/GUI args to settings dict:
    settings.update(params)

    settings.update(
        {
            "software_version": ver,
        }
    )

    # Write to terminal
    print("\n".join(f"{key}: {value}" for key, value in settings.items()) + "\n")

    """
    2. Format reference data and validate runinfo before counting
    """
    runinfo, p5_refdict, bc_refdict, p7_refdict = load_parse_reference_data(
        runinfo_path=settings["runinfo_path"],
        primer_path=settings["primer_path"],
        barcode_path=settings["barcode_path"],
        spike_path=settings.get("spike_path")
        if settings.get("append_spike_ref", False)
        else None,
        p7_path=settings["p7_path"] if settings["dualindex"] else None,
    )

    """
    3. Stream fastq to count barcode reads per index

    - Return {idx:{bc:count}}
    - See bcparse/ingest/fastq.py
    - mask returns N for base q<mask_qual
    - reads drop when mean phred<mean_qual
    - mismatches = # allowed errors
    """

    # Check for duplicate sample names before counting
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

        res_pfq = stream_to_countdict(
            settings=settings,
        )

        print("\n")

    """
    4. Main analysis: build one SeqRun from {idx_seq:{bc_seq:count}}
    """

    # Build parse-mode helper, then take the core SeqRun from it.
    parse_run = CountDictBuilder(
        countdict=res_pfq,
        runinfo=runinfo,
        bc_refdict=bc_refdict,
        p5_refdict=p5_refdict,
        p7_refdict=p7_refdict,
        settings=settings,
    ).build()

    seq_run = parse_run.seq_run

    print("\n")

    """
    5. Main output.
    """
    write_parse_workbook(builder=parse_run, settings=settings, runinfo_df=None)

    # In case of -i run, to examine:
    return res_pfq, seq_run, settings


# --- COMPILE MODE ---
def run_compile_mode(args):

    # Mode-specific imports:
    from bcparse.emit.workbooks import write_compile_workbook
    from bcparse.ingest.compile import build_runseries, read_base_csv
    from bcparse.ingest.xlsx import (
        AnalysisWorkbookParser,
        WorkbookParseError,
        XlsxPathManager,
    )

    # Define expected parameters
    # If running with --gui...
    if args.gui:
        from bcparse.gui.compile import CompileArgSelectionGUI

        # set variables returned from tkinter object,
        params = CompileArgSelectionGUI().run()
        print("GUI params:", params)

    # otherwise use CLI args.
    else:
        params = CompileSettings(
            xlsx_path=getattr(args, "xlsx_path"),
            out_path=getattr(args, "out_path"),
            base_csv_path=getattr(args, "base_csv_path"),
            out_prefix=getattr(args, "out_prefix"),
            ldist_samp_lvl=getattr(args, "ldist_samp_lvl"),
            ldist_group_lvl=getattr(args, "ldist_group_lvl"),
            contam_check=getattr(args, "contam_check"),
            dist_threshold=getattr(args, "dist_threshold"),
            collapse_to_parent=getattr(args, "collapse_to_parent"),
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
    settings["core_bc"] = (
        None  # For now, until I figure out how to detect stocks and do core search...
    )

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
    base_df = (
        read_base_csv(settings["base_csv_path"]) if settings["base_csv_path"] else None
    )
    if base_df is not None:
        print(
            f"Loaded base csv: {Path(settings['base_csv_path']).name} with {len(base_df)} row(s)."
        )

    print()

    # Iterate paths and obtain list of parsed xlsx data
    res_px = []

    for i, f in enumerate(xpm.files, start=1):
        print(f"Parsing workbook {i}/{len(xpm.files)}: {f.name}")
        try:
            px = AnalysisWorkbookParser(filepath=f)
            if not px.data_df.empty:
                res_px.append(px.analysis)
        # Exit on any parsing error, with informative message about which file caused the issue.
        except WorkbookParseError as e:
            raise SystemExit(f"Compile aborted while parsing workbook:\n{e}") from e

    print()

    # Compile set of parsed data to final output
    print("Building compiled series...")
    run_series = build_runseries(
        parsed_files=res_px,
        base_df=base_df,
        settings=settings,
    )

    print()
    print("Optional QC and output writing...")
    write_compile_workbook(run_series=run_series, settings=settings)
    print("Compile complete.")


# %% Main process function


def main():
    """
    Entry point for `python -m bcparse`.

    Parses CLI args and runs either Parse or Compile mode.
    Returns key objects so they can be accessed interactively
    (e.g., with `python -i -m bcparse ...`).
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


def cli() -> None:
    """
    Console-script entry point.
    """
    main()


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

# %%
