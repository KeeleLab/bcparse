# config.py

"""
Name:      config.py
Author:     CAG
Version:    1.7
Date:       2026/3/20

Description:
    Configure stock-specific settings for fastq read parsing.
    Values are added to settings dict in parse mode, and determine read handling logic in fastq_parser.py.
"""

# %% Imports
from pathlib import Path

# %% User-defined settings

# Reference data is packaged inside bcparse/ref.
package_dir = Path(__file__).resolve().parent
ref_data_dir = package_dir / "ref"


def ref_path(filename: str) -> str:
    path = ref_data_dir / filename
    if not path.is_file():
        raise FileNotFoundError(f"Missing packaged reference file: {path}")
    return str(path)


# Define stock settings
"""
General read structures -

M/VPX; 239M/M2/dual-index, NIRM, OptM, dGYM, SHIV:

    < 4bp tag :: (`tlen_p5`)bp idx :: `ref_p5`+ :: revcomp `ref_bc` :: revcomp (`tlen_bc`)bp barcode  >

    In 239M/M2, the biological data is reverse-complemented, but the index isn't, so the algotithm flow is:
    1. Identify `ref_p5` on raw sense strand
    2. Pull `tlen_p5` bases `upstream` of that
    3. `rdir` is rev, so reverse complement the sequence
    4. Identify `ref_bc` on revcomped new sense strand
    5. Pull `tlen_bc` bases `upstream` from the revcomped new sense strand

    In general, the index will always be apparent on the 5` end of the raw read.
    Revcomp occurs after identifying the index, when necessary for pulling region of interest.

    When running dualindex - modified read structure is as follows, prior to revcomp:
    < 239M/M2 read as above > :: < revcomp `ref_p7` :: revcomp (`tlen_p7`)bp p7 >

    NOTE that to match short ref barcodes from 239M/M2 p5-only characterization, provide
    a 'vpx_end' sequence which automatically pads barcodes shorter than 34bp as

    < 34 - length(bc_seq) bases 5' vpx_end 3' :: 5' bc_seq 3' >

    ... see parse_mode.py/parse_fastq_rec/match_extract_dual_targ().
    If vpx_end is not provided, the actual short sequence will be returned.

X/INT:

    < 4bp tag :: (`tlen_p5`)bp idx :: `ref_p5` = `ref_bc` :: (`tlen_bc`)bp barcode ... >

    More straightforward. Read is fwd-oriented, and the p5 index and barcode flank INT.

    In dual-index mode, p7 ends up revcomp - currently requires manual pre-processing via R1/R2 merge:
    Eg.
    ```
    vsearch \
        --fastq_mergepairs AVP082-SIVint-020626_S1_L001_R1_001.fastq.gz \
        --reverse AVP082-SIVint-020626_S1_L001_R2_001.fastq.gz \
        -- fastq_minovlen 10\
        --fastq_maxdiffs 10 \
        --fastq_allowmergestagger \
        --fastqout merged.test.fastq \
        --fastqout_notmerged_fwd /dev/null \
        --fastqout_notmerged_rev /dev/null
    ```

SL8_1-8:

    This is currently rarely used - built to rerun CTL_escapes.
    SL8 1-8 are disctinct from SL89-40, and I don't know why.

    Requires manual joining of R1/R2 via:
    ```
    vsearch \
        --fastq_join "$r1" \
        --reverse "$r2" \
        --join_padgap NNNN \
        --join_padgapq IIII \
        --fastqout - \
        --quiet |
        gzip -c >"${outdir}/${sample}.joined.fastq.gz"
    ```

    Read structure & treatment comparable to M/VPX, but with different p5/p7 ref targets & primer sets.

    < 4bp tag :: (`tlen_p5`)bp idx :: `ref_p5`+ :: revcomp `ref_bc` :: revcomp (`tlen_bc`)bp barcode :: revcomp `ref_p7` :: revcomp (`tlen_p7`)bp p7 >

    Primers added to P5/P7_primers.csv, require prefix "SL8" in runinfo, Eg. SL8.P5_1/SL8.P7_1.

"""

# Base read processing params
base_opts = {
    "M": {
        "rdir": "rev",
        "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",  # Start of VPR
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": ref_path("P5_primers.csv"),
        "append_spike_ref": False,
        "dualindex": False,
        "use_core_fallback": True,
        "core_bc": (
            12,
            10,
        ),  # Trim 12 nt from each end of the 34-mer; stored as 0-based start, length
    },
    "X": {
        "rdir": "fwd",
        "ref_p5": "TAAAAATTTTCGGGTCTATTAC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "TAAAAATTTTCGGGTCTATTAC",
        "tdir_bc": "downstream",
        "tlen_bc": 45,
        "primer_path": ref_path("P5_primers.csv"),
        "append_spike_ref": False,
        "dualindex": False,
        "core_bc": None,
    },
    "SL8_1-8": {
        "rdir": "rev",
        "ref_p5": "AGCTGAGAGAGGATTTCCTCCC",
        "tdir_p5": "upstream",
        "tlen_p5": 8,
        "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",
        "tdir_bc": "upstream",
        "tlen_bc": 34,
        "primer_path": ref_path("P5_primers.csv"),
        "append_spike_ref": False,
        "dualindex": False,
        "core_bc": None,
    },
}

# dual index settings
dual_index_opts = {
    "M": {
        "dualindex": True,
        "ref_p7": "CCTCCCCCTCCAGGACTAGCATAA",  # trunc end of VPX to p7 index
        "tdir_p7": "upstream",
        "tlen_p7": 8,
        "p7_path": ref_path("P7_primers.csv"),
        "fill_seq": "ACCTCCTCCTCCTCCCCCTCCAGGACTAGCATAA",  # VPX end to repair shorts
    },
    "X": {
        "dualindex": True,
        "ref_p7": "TGGATAGCAGTTCCCACATGGAGG",  # trunc end of INT to p7 index
        "tdir_p7": "downstream",
        "tlen_p7": 8,
        "p7_path": ref_path("P7_primers.csv"),
        "fill_seq": None,  # No fill for X series
    },
    # "SL8_1-8":{},
    # dual index is not currently applicable:
    # barcode = seq between ref_p7 and ref_bc, not fixed 34 nt upstream of ref_bc.
    # For joined SL8, this returns long construct sequence; the real 239M barcode would be the final 34 nt.
    # Using a longer ref_P7 would hurt due to potential mismatches.
    # Not sure of ideal fix; ignore for now. Single-index only.
}

controls_opts = {
    "Spike": {
        "spike_path": ref_path("Spike_reference.fasta"),
        "append_spike_ref": True,
    },
}

settings_opts = {
    "239M": {
        "base_profile": "M",
        "barcode_path": ref_path("239M_reference.fasta"),
    },
    "239M2": {
        "base_profile": "M",
        "barcode_path": ref_path("239M2_reference.fasta"),
    },
    "239X/INT": {
        "base_profile": "X",
        "barcode_path": ref_path("239X_reference.fasta"),
    },
    "M+M2": {
        "base_profile": "M",
        "barcode_path": ref_path("M_plus_M2_reference.fasta"),
    },
    "OptM": {
        "base_profile": "M",
        "barcode_path": ref_path("OptM_reference.fasta"),
    },
    "NIRM": {
        "base_profile": "M",
        "barcode_path": ref_path("NIRM_reference.fasta"),
    },
    "dGYM": {
        "base_profile": "M",
        "barcode_path": ref_path("dGYM_reference.fasta"),
    },
    "V67M": {
        "base_profile": "M",
        "barcode_path": ref_path("V67M_reference.fasta"),
    },
    "SHIV_ADE08M": {
        "base_profile": "M",
        "barcode_path": ref_path("SHIV_ADE08M_reference.fasta"),
    },
    "SHIV_174M": {
        "base_profile": "M",
        "barcode_path": ref_path("SHIV_174M_reference.fasta"),
    },
    "SHIV_224M": {
        "base_profile": "M",
        "barcode_path": ref_path("SHIV_224M_reference.fasta"),
    },
    "SHIV_304M": {
        "base_profile": "M",
        "barcode_path": ref_path("SHIV_304M_reference.fasta"),
    },
    "SHIV_1051M": {
        "base_profile": "M",
        "barcode_path": ref_path("SHIV_1051M_reference.fasta"),
    },
    "SHIV_1054M": {
        "base_profile": "M",
        "barcode_path": ref_path("SHIV_1054M_reference.fasta"),
    },
    "SL8_1-8_239M": {
        "base_profile": "SL8_1-8",
        "barcode_path": ref_path("239M_reference.fasta"),
    },
}

# List of available stocks
stocks = list(settings_opts.keys()) + ["Manual CLI"]

# %% Functions


# Define stock dict
def get_stock_settings(stock: str, dualindex: bool, append_spike_ref: bool = False):

    if stock in settings_opts:
        # Pull stock settings from dict
        settings_dict = settings_opts[stock].copy()
        base_profile = settings_dict.get("base_profile")

        # Update with base profile
        if base_profile not in base_opts:
            raise KeyError(f"Unknown base profile '{base_profile}' for stock {stock!r}")
        settings_dict.update(base_opts[base_profile].copy())

        if dualindex:
            # Update with di profile
            if base_profile not in dual_index_opts:
                raise KeyError(
                    f"Missing dual index profile for {stock!r}; run single index or add profile to dual_index_opts."
                )
            settings_dict.update(dual_index_opts[base_profile].copy())

        if append_spike_ref:
            if "Spike" not in controls_opts:
                raise KeyError("Unknown control option 'Spike'")
            settings_dict.update(controls_opts["Spike"].copy())

    else:
        print("Manual CLI mode — please enter all required parameters.")
        settings_dict = {
            "rdir": input("Enter read direction: "),
            "ref_p5": input("Enter p5 reference sequence: "),
            "tdir_p5": input("Enter p5 target direction: "),
            "tlen_p5": int(input("Enter p5 target length: ")),
            "ref_bc": input("Enter barcode reference sequence: "),
            "tdir_bc": input("Enter barcode target direction: "),
            "tlen_bc": int(input("Enter barcode target length: ")),
            "primer_path": input("Path/to/sequencing_primers.csv"),
            "barcode_path": input("Path/to/barcode_reference.fasta"),
            "spike_path": input(
                "Optional path/to/Spike_reference.fa (leave blank to skip): "
            ),
            "append_spike_ref": input("Append spike reference fasta? [y/n]: ")
            .lower()
            .startswith("y"),
            "dualindex": input("Run dual index mode? [y/n]: ").lower().startswith("y"),
            "di_profile": input("Profile for dual-index handling? [M/X]: "),
            "ref_p7": input("Enter p7 reference sequence: "),
            "tdir_p7": input("Enter p7 target direction: "),
            "tlen_p7": int(input("Enter p7 target length: ")),
            "p7_path": input("Path/to/p7_primers.csv"),
            "fill_seq": input(
                "Pad sequence for dual-index 3' repair - ask Charlie or `None`"
            ),
        }

        if not settings_dict["append_spike_ref"]:
            settings_dict["spike_path"] = None

    return settings_dict


# %% Versions
"""
1.7 - Added core_bc to base profile
1.6 - Dual index profile and config for 239X/INT, reduce redundant lines w/ base dicts.
1.5 - M/M2+Spike (Spike = SHIV174M renamed)
1.4 - Streamlined dualindex by adding seperate dict instead of redundant lines in all Ms.
To 1.3 - Various addition of new stocks
"""
