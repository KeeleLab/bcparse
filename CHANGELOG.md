# CHANGELOG

### [4.4.0] - Barcode discovery profiles
1. Added `M_discover` and `X_discover` parse profiles that extract and quantify
   barcodes without loading a stock barcode reference FASTA.
    - Name unmatched exact sequences as deterministic, run-wide uppercase `BC.N`
    - Apply the existing per-sample distance QC to all discovered barcodes
    - Emit `<run_name>_Discovery.xlsx`
    - Compile mode discovers both Analysis and Discovery workbooks
2. Promote RunInfo.normalize_input() to public, apply to emitted workbooks so "," in input values doesn't error out.

### [4.3.1] - New configs
1. Added TAT_SL8 epitope config, previous SL8 configuration reorganized in config.py
2. Modified parent-distance reporting to include `indel` in addition to substitutions

### [4.3.0] - tweaks
1. Containers;
  - comments and formatting, behaviour change - strict contracts,
    pipelines now own construction. Containers own validation.
    - runinfo.py v2.0.0 -> 2.0.1
    - runseries.py v2.0.0 -> 2.1.0
    - seqrun.py v1.0.0 -> v1.1.0
    - seqsamp.py v1.0.0 -> v1.1.0

2. Pipeline orchestration;
  - linear mode-owned parse and compile pipelines,
    retired countdict.py, fastq.py, compile.py, and xlsx.py as standalone ingest modules.
    - parse_mode.py / compile_mode.py

3. views.py v1.0.0 -> v1.1.0
  - fill source & sample_id on compiled runinfo

4. workbooks.py v2.0.0 -> v2.1.0
  - accomodating views.py change

5. parents.py v3.0.0 -> v3.1.0
  - Facilitate parse_mode.py w/ added arg
  - soothing pyright


### [4.2.0] - package refactor, uv distribution, new identity model
Restructure bcParse into a proper installable Python package with focused
submodules, plus substantial improvements to data identity, install/run
ergonomics, and code organization.

Package structure
- Flat utils/ layout replaced with bcparse/ package organized by domain:
  - containers/    runinfo, runseries, seqrun, seqsamp data models
  - ingest/        parse_mode and compile_mode linear pipelines
  - emit/          views, workbooks output rendering
  - qc/            collapse, flags, identity, parents (split from qc_ops.py)
  - gui/           parse, compile (split from monolithic gui.py)
  - lib/           bk_tree
  - settings.py    ParseSettings/CompileSettings dataclasses extracted
- Entry point moves from `python bcParse.py` to `python -m bcparse` or
  the installed `bcparse` command.

Installation and distribution
- uv-managed project: pyproject.toml, uv.lock, .python-version (3.13).
- Installable via uv tool / pipx / pip with isolated environments.
- Console script entry point: `bcparse` available on PATH after install.
- Reference data (ref/) packaged inside bcparse/ so installs are self-contained.
- Editable install supported for development (`uv tool install --editable .`).
- Dev dependencies: basedpyright, ruff, pytest, pandas-stubs.

Identity model overhaul
- New composite sample_id (group::name::date::run_number) is now the
  canonical compile-mode identity, replacing earlier idx_name-based dedup.
- New run_id (run_number::filename) provides stable cross-run references.
- Identity helpers consolidated in qc/identity.py with explicit validation:
  missing required fields raise clear errors instead of producing ambiguous
  ids.
- Compile mode now deduplicates on sample_id with provenance tracking
  (csv vs analysis, mtime, source order); duplicates are reported explicitly
  rather than silently overwritten.

QC subsystem split
- qc_ops.py (774 lines) split into focused modules:
  - qc/identity.py    id generation, normalization, long-df contract
  - qc/flags.py       proportion, above_cutoff, ranking, short_bc, x_idx
  - qc/parents.py     putative parent detection via BK-tree
  - qc/collapse.py    ambig and parent collapse operations
- RunSeries compile pipeline collapsed from ~15 staticmethods into
  build_runseries + _dedup + _core_qc.
- runinfo.py rewritten as a leaner class with explicit from_table /
  from_long_df constructors and clearer column contracts.

Other changes
- New SL8_1-8 stock configuration for CTL escape reruns.
- Renamed parameter: `collapse` → `collapse_ambig` for clarity (breaking).
- Type-checker (basedpyright) compliance throughout.
- Linter (ruff) configuration added.
- Per-sample progress reporting improved during qc_ops.
- Various small bug fixes in compile dedup, sample identity resolution,
  and date normalization.

Breaking changes from v4.1
- Installation now via uv/pipx/pip rather than running bcParse.py directly.
- Module imports restructured: `from utils.X import Y` → `from bcparse.X import Y`.
- The `--collapse` CLI flag is now `--collapse_ambig`.
- Compile dedup is now sample_id-based; runs with identical
  group/name/date/run_number that were previously treated as distinct may
  now be deduplicated. Review compile outputs against v4.1 if you have
  ambiguous metadata.
- Reference data (ref/) is now located inside the installed bcparse
  package rather than at a sibling path.
...

### [4.1.0] - growing pains! Current stable
1. Include ldist flags in all.csv compile output
2. Clean output dates to be YYYY-MM-DD
3. Added logic + checkbox to include Spike barcodes in parse mode, via settings and config.py. See --append_spike_ref and Spike entry in controls_opts for details.
4. Spike hits disallowed from above-cutoff box in sidelong tables
5. Normalize samp_group names as strings to avoid issues with 'number' animal names.

### [4.0.0] - Major overhaul of entire codebase, moving to class-based containers for runinfo, seqsamp, seqrun, and runseries.
0. All files drastically affected.
1. Change x_group to x_idx, and make it numeric for easier thresholding/idx hopping checks
2. Pre-flight parse mode sample name collision check, with informative error message
3. Optional collapse to parent in both parse and compile modes, default True. See settings and qc_ops.py 2.2 for details.
4. Resetting my comment log - previous edits can be found in deprecated branch bcParse_v3.42

### [3.4.3]
1. Dual-index X/INT config (config.py 1.5->1.6, parse_fastq.py 1.4->1.5)
2. Changed parent assignment logic for uniques: named > unique to highest named) (qc_ops.py 2.1->2.2)
3. To qc_ops
    - added informative suffix to parent name (qc_ops.py 2.1 -> 2.2 -> 2.3)
    - moved functions from parsefastq to qc_ops
    - big: core_bc search for unresolved uniques in flag_putative_parents
4. Dropped matrix backend outright, no need to carry forward

---

## Archived pre-git history

The following versions predate the current git repository. Source is
preserved in `bcParse_history_YYYYMMDD.tar.gz`.

### [3.4.2]
1. Add M/M2+Spike configs & ref data

### [3.4.1]
1. Normalize column datatypes in set_compiler._combine_with_base() to harden matrix pivot

### [3.4.0]
1. For analysis matrix output, modified to exclude above-cutoff children (omit second "block" in sidelong tables)
2. Added option to include a previously-compiled base 'all.csv' in compile mode
    - Behavior is to drop overlapping identical sample/runs from base and default to the new analysis
3. Fixed logger setup to be global (only I care about this)

### [3.3.2]
1. Fixed a major bug in combine_known_idx:
    - dist 1 between expected indexes (e.g. P5.40, P5.60) caused failed assignment!

### [3.3.1]
1. All I did was change primer index names from "VPX.P5/7.#" to "P5/7.#"

### [3.3.0]
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

### [3.2.0]
1. Repaired broken above-cutoff proportion normalization (was missing samp_group grouping factor)
2. Repaired contam_check missing flag

### [3.1.0]
1. Combined common methods across parse_countdata and set_compiler to utils.qc_ops
2. Modified compiler output naming variables
3. Added optional qc steps for compiled data incl. contam_check/full_mat
4. Inf input treated as 10 for cutoff determination, non-numeric throws error
5. Added SHIV Ms to config

### [3.0.0]
1. ==MAJOR== Added compiler mode w/ modules parse_xlsx.py and xlsx_compiler.py
2. Changed naming conventions of earlier modules for clarity.

### [2.5.3]
1. In cases where no named barcodes above cutoff are found, instead use nominal 1/input (min theoretical named) for uniques
2. Conditionally add per-index 'notes' column to runinfo
    - "ambig_children collapsed"
    - "no named above_cutoff"
3. Drop non-collapsed ambig uniques from results

### [2.5.2]
1. Updated config to run all "M" stocks with dual-index
2. Changed fastq streaming behavior to report count rather than progress bar - saves a lot of time

### [2.5.1]
1. Moved unique/above cutoff to rowgroup 0 for excel sample sheets
2. Recalculate rowgroup 0 props on ouput sample sheets to total 1

### [2.5.0]
1. Added modifiable quality params to adjust minimum req'd phred scores
    - mean_qual, mask_qual, and mismatches, passed through argparse or GUI
2. Changed date formats
3. Sort matrix output by date
4. New default run includes fuzzy matching and lowQ base masking
5. Changed argument handling across accessory scripts, all were reversioned
6. Substantially modified GUI
7. (Non-python) ./_misc_/MacroTemplate.xlsm included to recalculate proportions post-run.

### [2.4.3]
1. added samp_date to matrix output
2. modified i/o error message in response to a single-vs-dual index issue w/ identical sample names.
3. modified to handle low input samples such that "empty" indexes with no barcodes don't break pipe.

### [2.4.2]
1. Add columns to output runinfo:
    - n_reads, n_bc, n_reads_ac, n_bc_named, count-to-input_ratio (see io.py/write_output())
2. Modify output matrix:
    - Add sequencing input row to header
    - move bc_name and bc_seq to left

### [2.4.1]
1. Made output compatible with R-based PCRU concatination app (samp_date format problem)
2. Added 239X/INT, NIRM, OptM, and dGY stocks to config
3. Strip invisible whitespace from runinfo
4. Error check for known indexes - suggest runinfo format issue
5. Stress-test vs. large NIRM run

### [2.4.0]
1. Much time spent on making pyright happy (explicit type declarations, etc) - no functional changes

### [2.3.0]
1. removed --dualindex as argument and defined as per-stock setting in config.py (incl updates to gui and args)
2. added 239M/M2_dualindex as stock settings to accomodate 1.
3. BIG: for 239M/M2_dualindex, I extract VPX to VPR as barcode, whatever length, and fill to tlen_bc with vpx_end
        - This is to match previous versions 'short' barcodes
        - I think we should be reporting the actual extracted barcode instead, but whatever!
4. fastq_parser.py is fundamentally different now, incl. new match_extract_dual_targ()

### [2.2.0]:
1. parse_fastq_rec returned as instance instead of dict repr
2. p5/p7_seq and p5/p7_qual are held as independent attrs during parse_fastq_rec init, joined (or not) in process_fastq_stream()
3. bases < q20 are optionally substituted with N, resulting in idx and bc_seqs showing N, which can be leveraged for finer-grain analyses
4. --mask, --filt_mat_ac, --legacy_format added as optional boolean params, default False
5. GUI modified to indicate new params

### **Earlier versions preserved but untracked**

---

# Development Ideas

## Report M/M2 VPX overlap

- Keep dual-index `fill_seq` behavior: rebuild short observations into the expected 34 bp reference shape before barcode assignment.
- At dataframe QC, replace boolean `short_bc` with the inferred VPX suffix/barcode prefix overlap and report the remaining non-VPX sequence.
- Parse mode uses its known profile; compile mode infers a profile per sample from informative named barcodes and propagates it to Unique rows.
- Run overlap annotation only for M/M2-derived profiles; use `NA` for non-M or unresolved samples.

### Potential implementation

Replace `flag_shorts()` with a dataframe-level annotation. Parse mode passes its
known profile; compile mode supplies the inferred `barcode_profile` column.

```python
def flag_shorts(
    df: pd.DataFrame,
    *,
    fill_seq: str,
    base_profile: str | None = None,
    profile_col: str = "barcode_profile",
    min_overlap: int = 10,
) -> pd.DataFrame:
    df = df.copy()
    profiles = (
        pd.Series(base_profile, index=df.index, dtype="string")
        if base_profile is not None
        else df[profile_col].astype("string")
    )
    applies = profiles.isin({"M", "SL8_bc"})

    overlap = pd.Series(pd.NA, index=df.index, dtype="Int64")
    non_vpx = pd.Series(pd.NA, index=df.index, dtype="string")

    for idx, seq in df.loc[applies, "bc_seq"].items():
        seq = str(seq)
        matches = range(min(len(fill_seq), len(seq)), 0, -1)
        n_overlap = next(
            (n for n in matches if fill_seq[-n:] == seq[:n]),
            0,
        )
        n_overlap = n_overlap if n_overlap >= min_overlap else 0
        overlap.at[idx] = n_overlap
        non_vpx.at[idx] = seq[n_overlap:]

    df["short_bc"] = overlap
    df["non_vpx_seq"] = non_vpx
    return df
```

Compile ingest can resolve named barcodes against configured reference FASTAs,
then propagate the resulting barcode profile across each sample. Normalize
`SL8_bc` to `M` because it uses the 239M barcode reference biology.

```python
BARCODE_PROFILE = {
    "M": "M",
    "SL8_bc": "M",
    "X": "X",
    "SL8_epi": "SL8_epi",
}


def infer_barcode_profiles(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    wanted = set(df.loc[~df["bc_name"].str.startswith("Unique", na=False), "bc_name"])
    hits: dict[str, set[str]] = defaultdict(set)

    for stock in settings_opts.values():
        profile = BARCODE_PROFILE[stock["base_profile"]]
        with open(stock["barcode_path"], "rt") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    name = line[1:].strip()
                    if name in wanted:
                        hits[name].add(profile)

    name_profile = {
        name: next(iter(profiles))
        for name, profiles in hits.items()
        if len(profiles) == 1
    }

    def sample_profile(names: pd.Series):
        profiles = {name_profile[name] for name in names if name in name_profile}
        return next(iter(profiles)) if len(profiles) == 1 else pd.NA

    by_sample = df.groupby("sample_id")["bc_name"].agg(sample_profile)
    df["barcode_profile"] = df["sample_id"].map(by_sample).astype("string")
    return df
```

- Characterization mode for new stocks
  - Or at least: ability to run without barcode reference - for same purpose
- Coming quickly - ability to deal with degenerate primers (specifically HIVB...)
- Streamline running multiple parse runs - write arg to run a tsv of runinfos/fastqs/settings.
- Brandon wants a QC app... mode?
        - Develop this software into a more comprehensive analysis suite?
        - "live" mode, build a list of `bcsamp` objects, standardized figs, etc...

## Runinfo: sample-level passthrough column support to `RunInfo.from_long_df()`

`RunInfo.from_table()` already has a passthrough-column concept for extra
columns in the original runinfo sheet. `RunInfo.from_long_df()` currently builds
a minimal synthetic runinfo table from compiled/base long-format data, but it
does not preserve extra sample-level metadata columns that may have travelled
through the long dataframe.

Basic idea:

- Do not blindly copy every unknown long-format column into runinfo. Many extra
  long-format columns are barcode/count/QC fields, not sample metadata.
- Define a set of known long-format/internal columns to exclude, including run
  identifiers, sample identity columns, barcode fields, count/proportion fields,
  cutoff flags, ldist fields, contamination/QC annotations, and other derived
  analysis columns.
- Treat remaining columns as candidate sample-level passthrough columns.
- For each `sample_id` group, copy each candidate column only if it resolves to
  a single non-conflicting value within that sample, using the same
  `_single_value(..., sample_id=sample_id)` conflict check pattern already used
  for `idx_name`, `samp_group`, `samp_name`, `samp_date`, and `input`.
- Append accepted passthrough columns to the synthetic `raw_df` records and to
  `sample_rows`, after the required runinfo-derived columns.
- Keep synthetic bookkeeping columns such as `full_idx`, `notes`, and
  `sample_id` explicit so they do not collide with true user passthrough
  columns.

This would make reconstructed compile/base-CSV runinfo richer without changing
the core `from_long_df()` contract: input should still be one run-scoped
long-format dataframe, and passthrough values should still be sample-level
metadata, not per-barcode observations.

```python
sample_row_cols = [
    "idx_name",
    "samp_group",
    "samp_name",
    "samp_date",
    "input",
    "p7_name",
    "notes",
    "sample_id",
]

long_data_cols = {
    "run_id",
    "run_number",
    "run_name",
    "run_date",
    "filename",
    "sample_id",
    "idx_name",
    "samp_group",
    "samp_name",
    "samp_date",
    "input",
    "bc_name",
    "bc_seq",
    "bc_count",
    "proportion",
    "above_cutoff",
    "x_idx",
    "short_bc",
    "putative_parent",
    "ldist_samp_lvl",
    "ldist_group_lvl",
    "ldist_all",
    "multi_idx",
}

passthrough_cols = [c for c in df.columns if c not in long_data_cols]

# Inside the per-sample loop:
passthrough_values = {
    col: cls._single_value(sample_df[col], col, sample_id=sample_id)
    for col in passthrough_cols
}

sample_records.append(
    {
        "Run Number": run_number,
        "Animal": samp_group,
        "Sample": samp_name,
        "Date": samp_date,
        "Barcodes": barcode,
        "Input TOTAL PER BARCODE": input_val,
        "Input TOTAL PER WELL": "",
        "(F Barcode)": "" if p7_name is None else p7_name,
        "(cDNA)": "",
        "full_idx": idx_name,
        "notes": "from base_csv",
        "sample_id": sample_id,
        **passthrough_values,
    }
)

sample_rows = sample_rows[
    [c for c in sample_row_cols + passthrough_cols if c in sample_rows.columns]
].reset_index(drop=True)

return cls(
    ...,
    passthrough_cols=["full_idx", "notes", "sample_id"] + passthrough_cols,
    ...,
)
```

## Emit Reusable Run Settings

### Goal

Make frequent GUI-driven reruns easier by saving the exact final settings used for each analysis, then allowing a later run to reuse those settings.
This should work for both parse and compile mode, but parse mode is the higher-value first target because GUI parse reruns are common.

### Proposed UX

Every successful analysis writes a machine-readable settings sidecar to the output directory, for example:

```text
<run_name>_bcparse_settings.json
```

Later, the same settings can be reused from CLI:

```bash
bcparse --prev_settings path/to/run_bcparse_settings.json
```

Eventually, the GUI could add a "Load previous settings" button that populates the form from the same JSON file.

### Core Idea

Do not treat CLI arguments as the source of truth. Treat the final validated settings object as the source of truth.
Both GUI and CLI should produce the same mode-specific user settings object:

```text
GUI selections or CLI args
        ↓
ParseSettings / CompileSettings
        ↓
validate required user inputs
        ↓
resolve stock/config settings when relevant
        ↓
final settings dict/object used by the pipeline
        ↓
emit settings JSON to output directory
```

For parse mode, the final settings are the combination of:

- user selections: paths, stock, dual-index mode, QC flags, output options
- stock config: read direction, target references, target lengths, reference file paths
- run metadata: software version, settings schema version, creation timestamp

### Suggested JSON Shape

```json
{
  "bcparse_settings_version": 1,
  "mode": "parse",
  "software_version": "4.2 - 2026.05.27",
  "created_at": "2026-06-22 13:45:00",
  "settings": {
    "sample_path": "/path/to/sample.fastq.gz",
    "runinfo_path": "/path/to/runinfo.xlsx",
    "out_path": "/path/to/output",
    "stock": "239M",
    "dualindex": true,
    "mean_qual": 30,
    "mismatches": 1,
    "mask": false,
    "mask_qual": 20,
    "collapse_ambig": false,
    "dist_threshold": 1,
    "filt_mat_ac": true,
    "legacy_format": false,
    "collapse_to_parent": true,
    "append_spike_ref": false,
    "rdir": "rev",
    "ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
    "tdir_p5": "upstream",
    "tlen_p5": 8,
    "ref_bc": "ATGGAAGAAAGACCTCCAGAAAATGAAG",
    "tdir_bc": "upstream",
    "tlen_bc": 34,
    "primer_path": "/path/to/bcparse/ref/P5_primers.csv",
    "barcode_path": "/path/to/bcparse/ref/239M_reference.fasta"
  }
}
```

### Implementation Notes

Add a small settings IO helper, probably in a new module such as:

```text
bcparse/settings_io.py
```

Possible helpers:

```python
def write_settings_file(settings: dict, *, mode: str, out_path: str, prefix: str) -> Path:
    ...

def read_settings_file(path: str | Path) -> dict:
    ...
```

Keep this helper boring:

- use `json`
- use `Path`
- convert paths and non-JSON objects to strings before writing
- indent output for readability
- include a schema/version integer for future compatibility

### Where to Emit

Parse mode:

- after final settings are assembled in `run_parse_mode()`
- before or after `write_parse_workbook()`
- filename could use `seq_run.run_name` after the run is built, or a generic name before analysis starts

Best practical first version:

```text
<seq_run.run_name>_bcparse_settings.json
```

Compile mode:

```text
<out_prefix>_bcparse_settings.json
```

### Loading Previous Settings

Add a shared CLI argument:

```text
--prev_settings path/to/bcparse_settings.json
```

In each mode:

1. If `--prev_settings` is supplied, load the JSON.
2. Confirm the JSON `mode` matches the requested mode.
3. Use the stored `settings` dict as the starting settings.
4. Apply explicit current CLI overrides if supported.
5. Validate through `ParseSettings` or `CompileSettings`.

### Precedence Decision

Recommended precedence:

```text
built-in defaults < previous settings file < explicit current CLI args
```

This gives both convenient exact reruns and small overrides:

```bash
bcparse --prev_settings old_settings.json
```

```bash
bcparse --prev_settings old_settings.json --out_path new_output
```

Important caveat: the current argparse setup uses defaults for many options, so the code cannot always distinguish "user explicitly supplied this" from "argparse filled in the default." To make override precedence clean, argparse defaults should eventually move toward `None` or `argparse.SUPPRESS`, with defaults coming from `ParseSettings` / `CompileSettings`.

For a first GUI-oriented version, this caveat is less important: just emit the final settings file first. Loading and override semantics can come later.

### GUI Path

The GUI does not need a special implementation for emitting settings. It already returns a dict of selected values. That dict can become `ParseSettings`, then the resolved final settings can be emitted just like CLI settings.

Future GUI addition:

- add "Load previous settings"
- file dialog selects a JSON settings file
- populate GUI variables from `settings`
- keep output directory editable so reruns do not accidentally overwrite older output

### Safety Considerations

- Prefer absolute paths in emitted settings so reruns work from any launch directory.
- Consider also storing original display paths if relative paths are useful for readability.
- Validate that referenced input files still exist before rerunning.
- Do not silently overwrite previous output directories unless that is already normal bcparse behavior.
- Include `software_version` and `bcparse_settings_version` so old settings files can be detected later.

### Minimal First Pass

1. Add settings JSON writing after settings resolution.
2. Emit parse settings and compile settings sidecars.
3. Do not add `--prev_settings` yet.
4. Manually inspect one emitted JSON from a GUI run.

Then:

1. Add `--prev_settings`.
2. Load settings JSON.
3. Validate mode and required fields.
4. Run from loaded settings.

This keeps the first implementation small and useful even before the reload feature exists.

## Support Degenerate P5/P7 Target References

### Goal

Allow `ref_p5` and `ref_p7` to contain multiple acceptable concrete sequences.

The extracted index sequences remain concrete, so the existing index-reference and assignment logic does not need to change.

### Proposed Config Contract

Accept either a single sequence:

```python
"ref_p5": "CCAGAACCTCCACTACCCATTCATCC",
```

or multiple acceptable sequences:

```python
"ref_p5": (
    "CCAGAACCTCCACTACCCATTCATCC",
    "CCAGAACCTCCACTACCCGTTCATCC",
),
```

The same contract should apply to `ref_p7` and, where appropriate, `ref_bc`.

### Runtime Normalization

Normalize all target references to a tuple:

```python
TargetRef = str | tuple[str, ...]


def normalize_refs(value: TargetRef) -> tuple[str, ...]:
    return (value,) if isinstance(value, str) else value
```

This preserves compatibility with existing stock configurations.

### Pattern Compilation

Escape each concrete reference, join the alternatives with `|`, and compile them using the existing mismatch allowance:

```python
refs = normalize_refs(settings["ref_p5"])
alternatives = "|".join(regex.escape(ref) for ref in refs)

pat_p5 = regex.compile(
    rf"(?:{alternatives}){{s<={settings['mismatches']}}}"
)
```

This pattern means:

> Match any configured reference sequence while allowing at most
> `mismatches` base substitutions.

The existing `mismatches` setting already controls `{s<=...}` matching for `ref_p5`, `ref_p7`, and `ref_bc`, so no new mismatch parameter is required.

### Validation

For each target-reference collection:

- require at least one sequence;
- require nonempty sequences;
- normalize sequences to uppercase;
- allow only `A`, `C`, `G`, and `T`;
- reject duplicate alternatives;
- require all alternatives to have the same length.

Equal-length alternatives ensure that adjacent index-extraction boundaries remain consistent.

### Unchanged Behavior

No changes should be required for:

- extracting concrete P5/P7 index sequences;
- loading `P5_primers.csv` or `P7_primers.csv`;
- mapping `idx_seq` to `idx_name`;
- merging runinfo metadata;
- downstream barcode QC;
- barcode parent-distance checks.

### Distance Settings

The current workflow uses two related distance mechanisms:

```text
mismatches
├── ref_p5/ref_p7/ref_bc matching
│   └── substitution distance: {s<=mismatches}
└── observed index assignment
    └── Levenshtein distance

dist_threshold
└── downstream barcode parent/QC comparisons
```

Supporting multiple target-reference sequences can continue using the existing `mismatches` behavior.

### Focused Verification

Add smoke checks covering:

1. The original singular-string configuration still works.
2. Each configured `ref_p5` alternative extracts the expected P5 index.
3. Each configured `ref_p7` alternative extracts the expected P7 index.
4. An alternative containing one sequencing mismatch is accepted when `mismatches=1`.
5. The same mismatch is rejected when `mismatches=0`.
6. Invalid, duplicate, empty, or unequal-length references fail validation.
7. Extracted concrete index sequences still map to the expected sample.


# Current Architecture Notes

Date: 2026-06-25

This replaces the pre-consolidation Codex assessment. Resolved items from that
assessment were removed rather than kept as stale guidance.

Resolved since the original assessment:
- Parse and compile ingest now live in linear mode modules:
  - `bcparse/ingest/parse_mode.py`
  - `bcparse/ingest/compile_mode.py`
- Former standalone ingest modules were retired:
  - `countdict.py`
  - `fastq.py`
  - `compile.py`
  - `xlsx.py`
- `SeqSamp.from_table()` was removed; parse/compile modes now own tabular
  construction and call `SeqSamp(...)` directly.
- `SeqRun.from_long_df()` and flexible sample coercion were removed; parse and
  compile modes now materialize `SeqRun` objects explicitly.
- `RunSeries` construction is compile-mode owned and accepts one internal shape.
- `parents.py` now uses required `dist_threshold` access and has clean focused
  ruff/basedpyright checks.
- Mermaid diagrams were updated for the current mode-first/container-light flow.

Still-live cleanup themes:

1. Settings/default ownership
   - CLI, GUI, settings dataclasses, and runtime dicts still need a final call
     on where defaults and validation live.
   - Preferred direction remains: keep `ParseSettings`/`CompileSettings` as the
     user-facing config contract, convert to runtime dicts only where practical.

2. Long-dataframe contract boundaries
   - `normalize_long_df()` is still a broad repair/finalization function.
   - Future cleanup could split permissive external normalization from stricter
     internal contract validation.

3. Emit/view complexity
   - `emit/views.py` and `emit/workbooks.py` remain the densest user-facing
     formatting layer.
   - Sidelong table construction, matrix pivoting, date display normalization,
     and optional QC output are the main cleanup targets.

4. Legacy compatibility decisions
   - Decide which legacy workbook/input/output shapes still matter.
   - If old formats stay, quarantine them behind clearly named legacy adapters
     instead of letting them dominate the current path.

5. Known bug / output hardening
   - Parse output can still fail when raw runinfo input values contain comma
     separators, e.g. `Input TOTAL PER BARCODE = "300,000"`.
   - Fix should reuse normalized input values or coerce raw runinfo output fields
     before calculating workbook ratios.

6. Test and harness coverage
   - Add small characterization fixtures before deleting more emit/workbook or
     legacy compatibility code.
   - Minimum useful fixtures:
     - one FASTQ + runinfo parse fixture
     - one current-format Analysis workbook compile fixture

7. Packaging/dev-material cleanup
   - Revisit whether `_misc_` files should be shipped as package data.
   - Runtime references/templates should stay packaged; dev notes, mermaids,
     notebooks, R scratch, and workspace files probably should not ship.

8. Entry-point ownership
   - `bcparse/__main__.py` still owns CLI construction, orchestration, logging,
     interactive return objects, and user-facing help text.
   - It is workable, but future cleanup could split parser construction from
     parse/compile orchestration once behavior is stable.
