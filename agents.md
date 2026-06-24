# Agent Guidance

This repo is in an active refactor. The user's priority is not maximum
abstraction; it is linear legibility, strict contracts, and code that feels
owned rather than over-flexible.

## Permission Rule

Do not edit files unless the user explicitly asks for changes.

When the user is asking design questions, linter questions, review questions,
or "talk first" questions, stay read-only. Inspect code as needed, explain the
options, and wait for permission before patching.

The user often has files open and may be editing at the same time. Avoid
surprise edits, broad formatting churn, and unrelated cleanup.

## Refactor Direction

The desired direction is mode-first pipeline structure:

- `parse_mode.py` should tell the parse pipeline story.
- `compile_mode.py` should tell the compile pipeline story.
- Containers should be strict data holders and constructors from canonical data.
- Shared QC/identity/output helpers should stay separate when they are genuinely
  reusable.

The old `old_parse_countdata.py` style is an important reference point. Its
strength is an obvious, detailed, top-to-bottom pipeline in `__init__`, where
each major step either calls an internal helper method or a shared imported
function.

Prefer this shape:

```text
mode class __init__
  1. load inputs
  2. normalize/reference/validate
  3. build canonical long dataframe
  4. run QC
  5. construct containers at the end
  6. hand off to emit code
```

Avoid this shape:

```text
many containers and helper modules all accepting many input forms,
silently filling gaps, or hiding pipeline decisions in flexible constructors
```

## Container Philosophy

`RunInfo` is the model for the desired feel:

- It owns its construction machinery.
- `from_table()` builds from original runinfo/Analysis workbook tables.
- `from_long_df()` reconstructs synthetic runinfo from normalized long data.
- It has clear reasons for each constructor.

The other containers should move toward the same clarity.

### SeqSamp

`SeqSamp` should represent one logical sample. It should validate that its
dataframe belongs to one `sample_id`.

### SeqRun

Preferred future direction:

- `SeqRun.from_long_df()` becomes the canonical construction path for both parse
  and compile.
- Parse mode should build its final normalized long dataframe, then call
  `SeqRun.from_long_df(df, runinfo=runinfo, source="parse")`.
- Compile mode already reaches `SeqRun.from_long_df()` through
  `RunSeries.from_long_df()`.
- Remove unused flexibility such as accepting arbitrary iterables of `SeqSamp`.
- Remove empty/default object paths unless there is a real pipeline case.

### RunSeries

`RunSeries` has already been moved toward the desired style:

- It is compile-mode only.
- It is built from a non-empty final compiled long dataframe.
- `RunSeries.from_long_df()` splits by `run_id` and builds one `SeqRun` per run.
- Empty data should fail upstream in compile ingest.

## Mode Pipeline Sketch

### Parse Mode

Desired linear story:

```text
FASTQ/countdict + RunInfo.from_table()
  -> load reference data
  -> stream FASTQ into countdict
  -> assign observed indexes to expected indexes
  -> build barcode/count long dataframe
  -> merge runinfo sample metadata
  -> run per-sample QC
  -> normalize final long dataframe
  -> SeqRun.from_long_df(df, runinfo=runinfo, source="parse")
  -> emit parse outputs
```

### Compile Mode

Desired linear story:

```text
Analysis workbooks and optional base CSV
  -> discover/parse Analysis workbooks
  -> read optional base CSV
  -> normalize compile frames
  -> attach compile source provenance
  -> deduplicate by sample_id
  -> run core compile QC
  -> RunSeries.from_long_df(df, runinfo_by_run_id=...)
  -> ensure missing runinfo for base-CSV-derived runs
  -> emit compile outputs
```

## Current Provenance Context

Compile ingest now owns source provenance via an internal
`__compile_source_note` column.

`emit/views.py` should use that source note to populate a `source` column in the
compiled runinfo output. Do not append source text to `notes`.

Private `__compile_*` columns should not leak into user-facing CSV/XLSX outputs.

## Migration Strategy

Do not big-bang rewrite unless explicitly requested.

Prefer small, verifiable steps:

1. Make contracts stricter in place.
2. Align parse and compile on canonical long-dataframe construction.
3. Move orchestration into `parse_mode.py` and `compile_mode.py`.
4. Keep old function names temporarily as bridges.
5. Update imports after behavior is verified.
6. Remove shims and obsolete modules only after the new mode files are stable.

Useful temporary bridges:

```python
build_runseries(...)
load_parse_reference_data(...)
build_seqrun_from_countdict(...)
```

## Verification

For focused changes, run basedpyright on touched modules and nearby consumers,
for example:

```bash
uv run basedpyright bcparse/containers/seqrun.py bcparse/ingest/countdict.py
```

For behavior changes, prefer small smoke checks using the sample data in
`/Users/goodmanca/Documents/Scripts/BarcodeParser/_sampledata`.

When reporting command results, summarize the useful output in chat. The user
does not see tool output directly.
