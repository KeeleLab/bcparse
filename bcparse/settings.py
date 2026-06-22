# settings.py
"""
Name:       settings.py
Author:     CAG
Version:    1.0 
Date:       2026/05/26

Description:
    This module defines dataclasses to hold settings for parse and compile modes. 
    Descriptions are redundant with CLI args, but this structure allows for easier management of settings and validation.
"""

# %% Imports

from dataclasses import asdict, dataclass

# %% Classes

@dataclass
class ParseSettings:
    
    # Required
    sample_path:    str = ""    # Path to run.fastq
    runinfo_path:   str = ""    # Path to runinfo.xlsx
    out_path:       str = ""    # Path to output directory
    # Optional
    stock:              str = "239M"    # Default stock, check config.py for valid options
    dualindex:          bool = False    # Whether to parse dual-indexed samples (i7 and i5) or single-indexed (i7 only)
    mean_qual:          int = 30        # Minimum mean quality score for a read to be included in the output
    mismatches:         int = 1         # Mismatches allowed vs. reference sequences in read processing
    mask:               bool = False    # Mask low-quality bases
    collapse_ambig:     bool = False    # Collapse single-N children to parent rows
    mask_qual:          int = 20        # Quality threshold for masking bases
    dist_threshold:     int = 1         # Threshold for distance-based filtering
    filt_mat_ac:        bool = False    # Return above-cutoff matrix
    legacy_format:      bool = False    # Return legacy sidelong table format
    collapse_to_parent: bool = True     # Collapse ldist children to parent rows
    append_spike_ref:   bool = False    # Add 'spike' to reference

    # Helper methods for manual CLI ingestion and validation
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

    # Convert to plain dict for downstream use
    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class CompileSettings:

    # Required
    xlsx_path:      str = ""    # Path to dir with Analysis.xlsx files
    out_path:       str = ""    # Path to output directory
    out_prefix:     str = ""    # Prefix for output files (e.g. 'run1' to produce run1_compiled.csv)
    # Optional
    base_csv_path:      str | None = None    # Optional path to base CSV for merging with compiled data. 
    ldist_samp_lvl:     bool = False    # Rerun samp-lvl ldist
    ldist_group_lvl:    bool = False    # Run group-lvl ldist
    contam_check:       bool = False    # Run series ldist-based contamination check
    dist_threshold:     int = 1         # ldist to flag
    collapse_to_parent: bool = True     # Collapse ldist children to parent rows

    # Helper methods for manual CLI ingestion and validation
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

    # Convert to plain dict for downstream use
    def to_dict(self) -> dict:
        return asdict(self)


#%% Versions
"""
v1.0: Initial version with ParseSettings and CompileSettings dataclasses, including validation and dict conversion methods.
"""
