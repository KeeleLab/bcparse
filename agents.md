# Agent Guidance
The user's priority is not maximum abstraction; it is linear legibility, strict contracts, and code that feels
owned rather than over-flexible.

## Permission Rule

Do not edit files unless the user explicitly asks for changes.

When the user is asking design questions, linter questions, review questions,
or "talk first" questions, stay read-only. Inspect code as needed, explain the
options, and wait for permission before patching.

The user often has files open and may be editing at the same time. Avoid
surprise edits, broad formatting churn, and unrelated cleanup.

## General pipeline flow:

### Parse Mode

Desired linear story:

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
  -> build SeqRun objects
  -> RunSeries(...)
  -> ensure missing runinfo for base-CSV-derived runs
  -> emit compile outputs
```

## Migration Strategy

Do not big-bang rewrite unless explicitly requested.

Prefer small, verifiable steps:

1. Make contracts stricter in place.
2. Update imports after behavior is verified.
3. Remove shims and obsolete modules only after the new mode files are stable.

## Verification

For focused changes, run basedpyright on touched modules and nearby consumers,
for example:

```bash
uv run basedpyright bcparse/containers/seqrun.py bcparse/ingest/parse_mode.py
```

For behavior changes, prefer small smoke checks using the sample data in
`/Users/goodmanca/Documents/Scripts/BarcodeParser/_sampledata`.

When reporting command results, summarize the useful output in chat. The user
does not see tool output directly.
