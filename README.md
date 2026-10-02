# bcparse

`bcparse` is developed and maintained for the [Keele Lab - RES, ACVP, FNLCR](https://github.com/KeeleLab).

It is specifically designed for analysis of the lab's proprietary barcoded SIV assays,
and it supports two workflows:

- **Parse mode**: FASTQ plus runinfo workbook -> per-run `Analysis.xlsx` or `Discovery.xlsx` and long-format CSV
- **Compile mode**: one or more parse workbooks, optionally plus an existing compiled CSV -> compiled CSV and review workbook

The package includes the reference FASTA and primer CSV files in `bcparse/ref/`, so installed runs are self-contained.

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Parse Mode](#parse-mode)
- [Compile Mode](#compile-mode)
- [Changelog](#changelog)
- [License](#license)

## Installation

`bcparse` requires Python `>=3.13`.

We highly recommend uv package management:

```bash
# Standard
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or via homebrew
brew install uv
```

From a local checkout:

```bash
cd path/to/bcparse

# From top level
uv tool install .

# Or...
pipx install .

# Or...
pip install .

# To update, reinstall with --force (uv, pipx) or --upgrade (pip).
```

Install directly from GitHub:

```bash
uv tool install git+https://github.com/chazgoo/bcparse.git
pipx install git+https://github.com/chazgoo/bcparse.git
pip install git+https://github.com/chazgoo/bcparse.git

# For a specific tag, branch, or commit:
uv tool install git+https://github.com/chazgoo/bcparse.git@v4.4.0
```

Once installed, `bcparse` can be run from any directory:

```bash
bcparse --help
```

The launch directory dictates where GUI file dialogs open by default.

## Quick Start

Launch the parse GUI:

```bash
bcparse --gui
```

Launch the compile GUI:

```bash
bcparse --gui --compile
```

Parse one sequencing run from the command line:

```bash
bcparse \
  --sample_path path/to/sample.fastq.gz \
  --runinfo_path path/to/runinfo.xlsx \
  --out_path path/to/output \
  --stock 239M \
  --dualindex
```

Compile parsed analyses:

```bash
bcparse --compile \
  --xlsx_path path/to/analysis_dir \
  --out_path path/to/output \
  --out_prefix my_project
```

## Parse Mode

Parse mode processes one FASTQ or FASTQ.GZ file against one runinfo workbook.

Pipeline:

```text
FASTQ + RunInfo.from_table()
  -> load reference data
  -> stream FASTQ into countdict
  -> assign observed indexes to expected indexes
  -> build barcode/count long dataframe
  -> merge runinfo sample metadata
  -> run per-sample QC
  -> normalize final long dataframe
  -> build SeqSamp objects
  -> SeqRun(...)
  -> emit parse outputs
```

Inputs:

- `sample_path`: input FASTQ or gzipped FASTQ
- `runinfo_path`: runinfo Excel workbook
- `out_path`: output directory
- `stock`: stock-specific reference profile
- optional QC and formatting flags

Outputs:

- `<run_name>_concat.csv`: normalized long-format output
- `<run_name>_Analysis.xlsx`: reference-backed runinfo, sample tables, group
  matrices, and analysis settings
- `<run_name>_Discovery.xlsx`: the same workbook structure for discovery profiles

### Args
| Parameter                 | Required? | Default          | Description                                                            |
| ------------------------- | --------- | ---------------- | ---------------------------------------------------------------------- |
| `--sample_path`           | Yes       |                  | Path to input FASTQ/FASTQ.GZ file.                                     |
| `--runinfo_path`          | Yes       |                  | Path to runinfo Excel workbook.                                        |
| `--out_path`              | Yes       |                  | Output directory.                                                      |
| `--stock`                 | Yes       | `239M`           | Stock/reference profile to use.                                        |
| `--dualindex`             | No        | `False`          | Run in dual-index mode.                                                |
| `--mean_qual`             | No        | `30`             | Minimum allowed mean Phred score per read/sequence extract.            |
| `--mismatches`            | No        | `1`              | Number of allowed target mismatches during FASTQ parsing.              |
| `--mask`                  | No        | `False`          | Mask low-quality bases.                                                |
| `--mask_qual`             | No        | `20`             | Quality threshold for base masking.                                    |
| `--collapse_ambig`        | No        | `False`          | Collapse single-N ambiguous reads.                                     |
| `--legacy_format`         | No        | `False`          | Return output in legacy workbook/table format.                         |
| `--filt_mat_ac`           | No        | `False`          | Filter matrix output to above-cutoff barcodes only.                    |
| `--collapse_to_parent`    | No        | `True`           | Collapse putative child barcodes to parent rows.                       |
| `--no-collapse_to_parent` | No        | `False` override | Disable parent-collapse behavior.                                      |
| `--append_spike_ref`      | No        | `False`          | Append packaged Spike reference entries to selected barcode reference. |
| `--dist_threshold`        | No        | `1`              | Distance cutoff for parent checks.                                     |
| `--gui`                   | No        | `False`          | Use GUI for parse-mode file/option selection                           |
### Supported Stocks

Configured stock profiles currently include:

```text
M_discover
X_discover
239M
239M2
239X/INT
M+M2
OptM
NIRM
dGYM
V67M
SHIV_ADE08M
SHIV_174M
SHIV_224M
SHIV_304M
SHIV_1051M
SHIV_1054M
SL8_239M
TAT_SL8
```

Reference data is prepackaged in this repo.

### Barcode Discovery

`M_discover` and `X_discover` use the corresponding M or X extraction profile
without loading a stock barcode reference FASTA. Exact barcode sequences are
ranked across the run as uppercase `BC.1`, `BC.2`, and so on. Count ties are
resolved by sequence, and the same sequence keeps the same `BC.N` name in every
sample from that run.

Discovery barcodes participate in the normal per-sample full-sequence distance
checks. `BC.N` numbering is local to one run; use `bc_seq`, not the rank label,
as the durable molecular identity across runs.

The existing `--append_spike_ref` option remains available. Exact matches to an
appended reference retain their configured names, while unmatched sequences use
`BC.N` names and do not include exact reference matches in their rank numbering.

Example:

```bash
bcparse \
  --sample_path path/to/new_stock.fastq.gz \
  --runinfo_path path/to/runinfo.xlsx \
  --out_path path/to/output \
  --stock M_discover
```

`Manual CLI` also asks whether barcode discovery should be used. In manual
discovery, the barcode-reference path is omitted while the custom extraction
targets and primer paths are still required.

### Runinfo Workbook

Parse mode expects a 'runinfo.xlsx' Excel file matching this legacy format.

| Run Number              | Animal | Sample              | Date       | Barcodes | Input TOTAL PER BARCODE | Input TOTAL PER WELL |     | (F Barcode) | (cDNA) | Etc |
| ----------------------- | ------ | ------------------- | ---------- | -------- | ----------------------- | -------------------- | --- | ----------- | ------ | --- |
| MS3663441-300V2         | ML63   | 2WPI.Plasma.cDNA    | 2023-07-03 | P5.1     | 20000                   | 10000                |     | P7.39       | SL8R   |     |
|                         | RM07   | 30.Plasma.cDNA      | 2024-09-18 | P5.2     | 20000                   | 10000                |     | P7.39       | SL8R   |     |
|                         | RM07   | 53.PBMC.VifAluF.DNA | 2024-09-18 | P5.3     | 20000                   | 10000                |     | P7.39       | NA     |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
| Run Name                |        |                     |            |          |                         |                      |     |             |        |     |
| P7newIndex2-test-011325 |        |                     |            |          |                         |                      |     |             |        |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
| Date                    |        |                     |            |          |                         |                      |     |             |        |     |
| 1/13/25                 |        |                     |            |          |                         |                      |     |             |        |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
| Regular Taq             |        |                     |            |          |                         |                      |     |             |        |     |
| 300-V2 Nano Kit         |        |                     |            |          |                         |                      |     |             |        |     |
| VPX                     |        |                     |            |          |                         |                      |     |             |        |     |
|                         |        |                     |            |          |                         |                      |     |             |        |     |
| One-way                 |        |                     |            |          |                         |                      |     |             |        |     |
This runinfo format is required.

- Run-level metadata is read from the first column.
    - `Run Number`, `Run Name`, and `Date` must be filled in using the row-wise layout shown in the example.
    - Any additional rows below those fields are optional. They are carried through to the output runinfo sheet but are not used by the parser.

- Sample identity requires `Animal`, `Sample`, and `Date`.
    - Each `Animal` / `Sample` / `Date` combination must be unique within a run.
    - Replicates should be given distinct sample names or numbers.

- `Barcodes` lists the P5 sequencing index assigned to each sample.
    - For dual-index runs, `(F Barcode)` supplies the P7 index.
    - In dual-index mode, `Barcodes` and `(F Barcode)` are combined to form the expected P5/P7 index identity.

- `Input TOTAL PER BARCODE` is required.
    - It is used to calculate barcode abundance cutoffs.

- `Input TOTAL PER WELL` and `(cDNA)` are accepted for compatibility but are not used by `bcparse`.

- Additional metadata columns may be added after the required columns.
    - These columns are passed through to the output runinfo sheet.

## Compile Mode

Compile mode combines parse-mode `*Analysis.xlsx` and `*Discovery.xlsx` files
into a cross-run dataset.

Pipeline:

```text
Analysis/Discovery workbooks and optional base CSV
  -> discover/parse workbooks
  -> read optional base CSV
  -> normalize compile frames
  -> attach compile source provenance
  -> deduplicate by sample_id
  -> run core compile QC
  -> build SeqRun objects
  -> RunSeries(...)
  -> ensure missing runinfo for base-CSV-derived runs
  -> emit compile outputs
```

Inputs:

- `xlsx_path`: directory containing `*Analysis.xlsx` and/or `*Discovery.xlsx` workbooks
- `out_path`: output directory
- `out_prefix`: compiled output name prefix
- optional `base_csv_path`
- optional compile-side QC flags

Outputs:

- `<out_prefix>_all.csv`: full compiled long-format table
- `<out_prefix>_above_cutoff.xlsx`: above-cutoff review workbook with runinfo,
  sample tables, matrices, and optional QC sheets

Compile mode deduplicates by `sample_id`. When the same sample is present in
multiple inputs, parsed workbooks are preferred over base CSV rows,
then newer workbook modification time and input order break remaining ties.

### Args
| Parameter                 | Required? | Default          | Description                                                                         |
| ------------------------- | --------- | ---------------- | ----------------------------------------------------------------------------------- |
| `--compile`               | Yes       | `False`          | Run compile mode instead of parse mode.                                             |
| `--xlsx_path`             | Yes       |                  | Path to directory containing `*Analysis.xlsx` and/or `*Discovery.xlsx` files.       |
| `--out_path`              | Yes       |                  | Output directory.                                                                   |
| `--out_prefix`            | Yes       |                  | Prefix for compiled output files.                                                   |
| `--base_csv_path`         | No        | `None`           | Optional path to a previous compiled `*_all.csv` to merge with new parsed analyses. |
| `--ldist_samp_lvl`        | No        | `False`          | Rerun sample-level distance checks.                                                 |
| `--ldist_group_lvl`       | No        | `False`          | Run group-level distance checks.                                                    |
| `--contam_check`          | No        | `False`          | Run contamination check and emit extra workbook sheets.                             |
| `--collapse_to_parent`    | No        | `True`           | Collapse putative child barcodes to parent rows.                                    |
| `--no-collapse_to_parent` | No        | `False` override | Disable parent-collapse behavior.                                                   |
| `--dist_threshold`        | No        | `1`              | Distance cutoff for parent checks.                                                  |
| `--gui`                   | No        | `False`          | Use GUI for compile-mode file/option selection.                                     |
|                           |           |                  |                                                                                     |
## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

Add a `LICENSE` file before publishing this repository.
