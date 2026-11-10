# screen_mode.py
"""
Name:       screen_mode.py
Author:     CAG
Version:    1.3.0
Date:       2026/11/10

Run-level FASTQ screening across configured read profiles and stocks.

This module is intentionally linear, ScreenMode.__init__ describes action.

The goal is to identify read-profile and stock composition before RunInfo,
sample assignment, per-sample QC, or SeqRun construction.

Major classes and helpers:

    ScreenMode
        Owns profile-parser construction, stock/control reference loading,
        FASTQ screening, optional barcode-distance matching, and output tables.

    ProfileHit
        Small immutable record for one quality-qualified profile extraction.

    fq_profile_parser / stream_fastq
        Shared parse-mode machinery for structural target extraction and FASTQ
        streaming.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.metadata import version as package_version
from pathlib import Path
from typing import cast

import pandas as pd
from tqdm import tqdm

from bcparse.config import (
    get_screen_profile_settings,
    get_screen_stock_settings,
)
from bcparse.modes.parse_mode import fq_parsed_rec, fq_profile_parser, stream_fastq
from bcparse.lib.bk_tree import BKTree, MinAllResult, build_ref_index
from bcparse.qc.parents import parent_suffix_sub
from bcparse.settings import ScreenSettings

# %% ScreenMode support

@dataclass(frozen=True)
class ProfileHit:
    """One quality-qualified profile, structure, and extracted barcode."""

    profile: str
    structure: str
    barcode: str


@dataclass(frozen=True)
class StockMatch:
    """Reference barcode names, stock owners, and distance for one match."""

    barcode_names: tuple[str, ...]
    ref_matches: tuple[str, ...]
    stocks: tuple[str, ...]
    distance: int | None


# %% ScreenMode

class ScreenMode:
    """
    Top-level screen-mode orchestration for one FASTQ.

    Call site:
        bcparse.__main__.run_screen_mode() builds this class as:
            screen_run = ScreenMode(settings=screen_settings)

    Inputs:
        ScreenSettings object from CLI/user configuration:
            - sample FASTQ path
            - output path
            - anchor mismatch and quality thresholds
            - optional barcode distance threshold

    Outputs:
        summary_df:
            Terminal read outcomes which sum to total FASTQ reads.
        profiles_df:
            Counts and proportions for uniquely classified read profiles.
        stocks_df:
            Counts and proportions for uniquely assigned stocks/controls.
        ambiguous_df:
            Profile or stock ties retained for explicit review.
        barcode_dfs:
            Parse-style long barcode tables keyed by classified profile.
        ambiguous_barcodes_df:
            Barcode-level detail for profile and stock ties.
    """

    settings: ScreenSettings  # Original user/CLI screen config.
    runtime_settings: dict[str, object]  # Emitted screen settings/version.

    profile_parsers: dict[str, dict[str, fq_profile_parser]]  # Profile -> structure -> parser.
    stock_reference_paths: dict[str, dict[str, str]]  # Profile -> stock/control -> FASTA path.
    stock_lookups: dict[str, dict[str, tuple[tuple[str, str], ...]]]  # Profile -> bc_seq -> (bc_name, stock).
    stock_trees: dict[str, BKTree]  # Optional per-profile barcode distance indexes.

    outcome_counts: Counter[str]  # Terminal read outcome counts.
    profile_counts: Counter[tuple[str, str]]  # (profile, structure) counts.
    extracted_barcode_counts: Counter[tuple[str, str, str]]  # (profile, structure, bc_seq) counts.
    ambiguous_profile_hit_counts: Counter[tuple[str, str, str, str]]  # (profiles, profile, structure, bc_seq) counts.
    stock_counts: Counter[tuple[str, str, int]]  # (profile, stock, ldist) counts.
    ambiguous_counts: Counter[tuple[str, str, int | None]]  # (kind, labels, ldist) counts.
    barcode_counts: Counter[tuple[str, str, str, tuple[str, ...], tuple[str, ...], tuple[str, ...], int | None]]  # Profile barcode rows.
    ambiguous_barcode_counts: Counter[tuple[str, str, str, str, str, tuple[str, ...], tuple[str, ...], tuple[str, ...], int | None]]  # Ambiguous barcode rows.

    summary_df: pd.DataFrame  # Terminal outcome summary table.
    profiles_df: pd.DataFrame  # Profile/structure composition table.
    stocks_df: pd.DataFrame  # Uniquely assigned stock/control table.
    ambiguous_df: pd.DataFrame  # Profile and stock tie table.
    barcode_dfs: dict[str, pd.DataFrame]  # Profile -> long barcode table.
    ambiguous_barcodes_df: pd.DataFrame  # Barcode-level profile/stock ties.

    def __init__(self, *, settings: ScreenSettings) -> None:

        # 1. Validate ScreenSettings and build emitted runtime settings.
        settings.validate()
        self.settings = settings
        self.runtime_settings = settings.to_dict()
        self.runtime_settings["software_version"] = package_version("bcparse")

        # 2. Build independent single/dual parsers for every read profile.
        self.profile_parsers = self._build_profile_parsers()

        # 3. Load stock/control references into per-profile sequence lookups.
        self.stock_reference_paths = get_screen_stock_settings()
        self.stock_lookups = self._load_stock_lookups()

        # 4. Optionally build per-profile barcode-distance indexes.
        self.stock_trees = self._build_stock_trees()

        # 5. Initialize read, extracted-barcode, and result counters.
        self.outcome_counts = Counter()
        self.profile_counts = Counter()
        self.extracted_barcode_counts = Counter()
        self.ambiguous_profile_hit_counts = Counter()
        self.stock_counts = Counter()
        self.ambiguous_counts = Counter()
        self.barcode_counts = Counter()
        self.ambiguous_barcode_counts = Counter()

        # 6. Stream and classify the FASTQ.
        self._screen_fastq()

        # 7. Resolve each distinct profile/barcode sequence once.
        self._resolve_barcodes()

        # 8. Build workbook-ready output tables.
        self._build_output_frames()

    def _build_profile_parsers(
        self,
    ) -> dict[str, dict[str, fq_profile_parser]]:
        """
        Build single/dual parser variants for every configured read profile.

        Screen-level quality and mismatch settings override the packaged profile
        defaults without changing config-owned structural settings.
        """
        # Settings shared across every structural profile.
        common = {
            "mismatches": self.settings.mismatches,
            "mask": self.settings.mask,
            "mask_qual": self.settings.mask_qual,
        }

        # Compile each profile/structure once for reuse across all FASTQ records.
        parsers: dict[str, dict[str, fq_profile_parser]] = {}
        for profile, variants in get_screen_profile_settings().items():
            parsers[profile] = {}
            for structure, profile_settings in variants.items():
                parser_settings = {**profile_settings, **common}
                parsers[profile][structure] = fq_profile_parser(parser_settings)
        return parsers

    @staticmethod
    def _read_fasta_records(path: str) -> dict[str, str]:
        """Read one FASTA into a barcode-name -> sequence mapping."""
        records: dict[str, str] = {}
        current_name: str | None = None
        current: list[str] = []

        # Join multiline FASTA records while retaining their barcode names.
        with open(path, "r") as fasta:
            for raw_line in fasta:
                line = raw_line.strip()
                if line.startswith(">"):
                    if current_name is not None:
                        records[current_name] = "".join(current)
                    current_name, current = line[1:], []
                elif line:
                    current.append(line)
        if current_name is not None:
            records[current_name] = "".join(current)
        return records

    def _load_stock_lookups(
        self,
    ) -> dict[str, dict[str, tuple[tuple[str, str], ...]]]:
        """
        Build profile -> sequence -> (barcode name, stock) lookups.

        Multiple records are retained when the same exact barcode belongs to
        more than one reference set or has more than one reference name.
        """
        lookups: dict[str, dict[str, tuple[tuple[str, str], ...]]] = {}
        for profile, stock_paths in self.stock_reference_paths.items():
            references: defaultdict[str, set[tuple[str, str]]] = defaultdict(set)

            # Collect every named stock reference for each exact sequence.
            for stock, path in stock_paths.items():
                for barcode_name, sequence in self._read_fasta_records(path).items():
                    references[sequence].add((barcode_name, stock))

            # Freeze and sort records for deterministic downstream reporting.
            lookups[profile] = {
                sequence: tuple(sorted(records))
                for sequence, records in references.items()
            }
        return lookups

    def _build_stock_trees(self) -> dict[str, BKTree]:
        """Build per-profile barcode-distance indexes only when requested."""
        # Exact-only screening needs no distance index.
        if self.settings.dist_threshold == 0:
            return {}

        # Index each unique reference sequence once per structural profile.
        trees: dict[str, BKTree] = {}
        for profile, lookup in self.stock_lookups.items():
            sequences = sorted(lookup)
            if sequences:
                trees[profile] = build_ref_index(sequences, sequences)
        return trees

    def _match_stocks(
        self,
        *,
        profile: str,
        barcode: str,
    ) -> StockMatch:
        """
        Resolve one extracted barcode to nearest named stock/control records.

        Policy:
        1) Exact reference matches always win at ldist 0.
        2) Otherwise retain all equal-nearest references within dist_threshold.
        3) Return no owners when no permitted match exists.

        _resolve_barcodes() calls this once per distinct profile/barcode pair.
        """
        # Prefer exact reference ownership before any BK-tree query.
        lookup = self.stock_lookups[profile]
        exact = lookup.get(barcode)
        if exact is not None:
            return StockMatch(
                barcode_names=tuple(sorted({name for name, _ in exact})),
                ref_matches=tuple(sorted({name for name, _ in exact})),
                stocks=tuple(sorted({stock for _, stock in exact})),
                distance=0,
            )

        # No tree means exact-only mode or no references for this profile.
        tree = self.stock_trees.get(profile)
        if tree is None:
            return StockMatch(
                barcode_names=tuple(),
                ref_matches=tuple(),
                stocks=tuple(),
                distance=None,
            )

        # Retain all sequences tied at the smallest permitted edit distance.
        nearest_sequences, distance = cast(
            MinAllResult,
            tree.nearest(
                barcode,
                max_radius=self.settings.dist_threshold,
                exclude_same_id=False,
                tie="all",
                scope="min",
            ),
        )
        # Preserve each tied reference sequence long enough to describe its edit.
        reference_records = {
            (name, stock, str(sequence))
            for sequence in nearest_sequences
            for name, stock in lookup[str(sequence)]
        }
        return StockMatch(
            barcode_names=tuple(
                sorted({name for name, _, _ in reference_records})
            ),
            ref_matches=tuple(
                sorted(
                    {
                        f"{name}_{parent_suffix_sub(reference, barcode)}"
                        for name, _, reference in reference_records
                    }
                )
            ),
            stocks=tuple(sorted({stock for _, stock, _ in reference_records})),
            distance=distance,
        )

    @staticmethod
    def _has_extract(result: fq_parsed_rec) -> bool:
        """Return whether both index and biological barcode were extracted."""
        return hasattr(result, "p5_seq") and hasattr(result, "bc_seq")

    def _parse_profile(
        self,
        *,
        profile: str,
        seq: str,
        qual: str,
    ) -> tuple[fq_parsed_rec, str]:
        """
        Parse one read against one profile, preferring a detected dual structure.

        Profiles without a detected P7 target fall back to their single-index
        parser variant.
        """
        variants = self.profile_parsers[profile]

        # Try dual first so joined reads use P7-aware barcode extraction/repair.
        dual_parser = variants.get("dual")
        if dual_parser is not None:
            dual_result = dual_parser.parse(seq=seq, qual=qual)
            if hasattr(dual_result, "p7_seq"):
                return dual_result, "dual"

        # No configured/detected P7: use the single-index structure.
        single_result = variants["single"].parse(seq=seq, qual=qual)
        return single_result, "single"

    def _screen_fastq(self) -> None:
        """
        Stream one FASTQ and count its quality-qualified profile extractions.

        Structural failures and profile ties receive terminal outcomes here.
        Unique-profile barcode sequences are retained for _resolve_barcodes().
        """
        for _, seq, _, qual in stream_fastq(self.settings.sample_path):
            self.outcome_counts["total_reads"] += 1
            extracted: list[tuple[fq_parsed_rec, str, str]] = []

            # Attempt structural extraction against every configured profile.
            for profile in self.profile_parsers:
                result, structure = self._parse_profile(
                    profile=profile,
                    seq=seq,
                    qual=qual,
                )
                if self._has_extract(result):
                    extracted.append((result, profile, structure))

            # No profile produced both an index and biological barcode.
            if not extracted:
                self.outcome_counts["profile_not_found"] += 1
                continue

            # Retain only extractions meeting P5 and barcode quality thresholds.
            hits = [
                ProfileHit(
                    profile=profile,
                    structure=structure,
                    barcode=result.bc_seq,
                )
                for result, profile, structure in extracted
                if result.p5_qual >= self.settings.mean_qual
                and result.bc_qual >= self.settings.mean_qual
            ]

            # Structural extraction succeeded, but all candidates failed quality.
            if not hits:
                self.outcome_counts["quality_filtered"] += 1
                continue

            # Do not double-count a read which qualifies for multiple profiles.
            if len(hits) > 1:
                profiles = tuple(sorted({hit.profile for hit in hits}))
                labels = "/".join(profiles)
                self.outcome_counts["profile_ambiguous"] += 1
                self.ambiguous_counts[("profile", labels, None)] += 1

                # Retain each candidate extraction because profiles may extract
                # different barcode sequences from the same ambiguous read.
                for hit in sorted(hits, key=lambda item: item.profile):
                    self.ambiguous_profile_hit_counts[
                        (
                            labels,
                            hit.profile,
                            hit.structure,
                            hit.barcode,
                        )
                    ] += 1
                continue

            # One profile qualified: retain its extracted sequence for resolution.
            hit = hits[0]
            self.profile_counts[(hit.profile, hit.structure)] += 1
            self.extracted_barcode_counts[
                (hit.profile, hit.structure, hit.barcode)
            ] += 1

    def _resolve_barcodes(self) -> None:
        """Match each distinct profile/barcode pair and weight by read count."""
        observed_barcodes = {
            (profile, barcode)
            for profile, _, barcode in self.extracted_barcode_counts
        }
        observed_barcodes.update(
            (profile, barcode)
            for _, profile, _, barcode in self.ambiguous_profile_hit_counts
        )

        # Distance work happens here exactly once per profile/barcode pair.
        matches: dict[tuple[str, str], StockMatch] = {}
        for profile, barcode in tqdm(
            sorted(observed_barcodes),
            desc="Matching barcode references",
            unit="bc",
            ascii="-=",
        ):
            matches[(profile, barcode)] = self._match_stocks(
                profile=profile,
                barcode=barcode,
            )

        # Resolve barcodes from reads with one qualified structural profile.
        for key, count in self.extracted_barcode_counts.items():
            profile, structure, barcode = key
            match = matches[(profile, barcode)]
            self.barcode_counts[
                (
                    profile,
                    structure,
                    barcode,
                    match.barcode_names,
                    match.ref_matches,
                    match.stocks,
                    match.distance,
                )
            ] += count

            # Keep unknowns and ties explicit instead of forcing assignment.
            if not match.stocks:
                self.outcome_counts["stock_unknown"] += count
            elif len(match.stocks) > 1:
                labels = "/".join(match.stocks)
                self.outcome_counts["stock_ambiguous"] += count
                self.ambiguous_counts[
                    ("stock", labels, match.distance)
                ] += count
                self.ambiguous_barcode_counts[
                    (
                        "stock",
                        profile,
                        profile,
                        structure,
                        barcode,
                        match.barcode_names,
                        match.ref_matches,
                        match.stocks,
                        match.distance,
                    )
                ] += count
            else:
                if match.distance is None:
                    raise RuntimeError(
                        "Assigned stock match is missing its distance."
                    )
                self.outcome_counts["stock_assigned"] += count
                stock = next(iter(match.stocks))
                self.stock_counts[(profile, stock, match.distance)] += count

        # Resolve each candidate extraction from profile-ambiguous reads.
        for key, count in self.ambiguous_profile_hit_counts.items():
            candidate_profiles, profile, structure, barcode = key
            match = matches[(profile, barcode)]
            self.ambiguous_barcode_counts[
                (
                    "profile",
                    candidate_profiles,
                    profile,
                    structure,
                    barcode,
                    match.barcode_names,
                    match.ref_matches,
                    match.stocks,
                    match.distance,
                )
            ] += count

    @staticmethod
    def _proportion(count: int, denominator: int) -> float:
        """Return a zero-safe proportion."""
        return 0.0 if denominator == 0 else count / denominator

    def _build_barcode_frames(self, *, total_reads: int) -> None:
        """Build parse-style per-profile and ambiguous barcode tables."""
        columns = [
            "bc_name",
            "bc_count",
            "proportion",
            "bc_seq",
            "ref_match",
            "stock",
            "ldist",
            "index_structure",
            "proportion_all_reads",
        ]
        self.barcode_dfs = {}

        # Build one long table for each confidently classified profile.
        profiles = sorted({key[0] for key in self.barcode_counts})
        for profile in profiles:
            rows = [
                {
                    "bc_name": (
                        barcode_names[0]
                        if distance == 0 and len(ref_matches) == 1
                        else "Multi"
                        if len(ref_matches) > 1
                        else "Unique"
                    ),
                    "bc_count": count,
                    "bc_seq": barcode,
                    "ref_match": (
                        ""
                        if distance == 0 and len(ref_matches) == 1
                        else "/".join(ref_matches)
                    ),
                    "stock": "/".join(stocks),
                    "ldist": distance,
                    "index_structure": structure,
                }
                for (
                    row_profile,
                    structure,
                    barcode,
                    barcode_names,
                    ref_matches,
                    stocks,
                    distance,
                ), count in self.barcode_counts.items()
                if row_profile == profile
            ]
            df = pd.DataFrame(rows)
            profile_total = int(df["bc_count"].sum())

            # Rank non-exact and multi-reference observations independently.
            for prefix in ("Unique", "Multi"):
                observed = df["bc_name"] == prefix
                sequence_counts = (
                    df.loc[observed]
                    .groupby("bc_seq")["bc_count"]
                    .sum()
                    .to_dict()
                )
                ranked_names = {
                    sequence: f"{prefix}.{rank}"
                    for rank, sequence in enumerate(
                        sorted(
                            sequence_counts,
                            key=lambda sequence: (
                                -sequence_counts[sequence],
                                sequence,
                            ),
                        ),
                        start=1,
                    )
                }
                df.loc[observed, "bc_name"] = df.loc[
                    observed, "bc_seq"
                ].map(ranked_names)
            df.loc[df["stock"] == "", "stock"] = "Unknown"

            df["proportion"] = df["bc_count"] / profile_total
            df["proportion_all_reads"] = (
                0.0 if total_reads == 0 else df["bc_count"] / total_reads
            )
            self.barcode_dfs[profile] = df.loc[:, columns].sort_values(
                ["bc_count", "bc_seq"],
                ascending=[False, True],
            ).reset_index(drop=True)

        ambiguous_columns = [
            "ambiguity_type",
            "profiles",
            "profile",
            "index_structure",
            "bc_name",
            "bc_count",
            "proportion",
            "bc_seq",
            "ref_match",
            "stocks",
            "ldist",
            "proportion_all_reads",
        ]
        ambiguous_base_columns = [
            column
            for column in ambiguous_columns
            if column not in {"proportion", "proportion_all_reads"}
        ]
        ambiguous_rows = [
            {
                "ambiguity_type": kind,
                "profiles": profiles,
                "profile": profile,
                "index_structure": structure,
                "bc_name": (
                    barcode_names[0]
                    if distance == 0 and len(ref_matches) == 1
                    else "Multi"
                    if len(ref_matches) > 1
                    else "Unique"
                ),
                "bc_count": count,
                "bc_seq": barcode,
                "ref_match": (
                    ""
                    if distance == 0 and len(ref_matches) == 1
                    else "/".join(ref_matches)
                ),
                "stocks": "/".join(stocks),
                "ldist": distance,
            }
            for (
                kind,
                profiles,
                profile,
                structure,
                barcode,
                barcode_names,
                ref_matches,
                stocks,
                distance,
            ), count in self.ambiguous_barcode_counts.items()
        ]
        self.ambiguous_barcodes_df = pd.DataFrame(
            ambiguous_rows,
            columns=ambiguous_base_columns,
        )
        if self.ambiguous_barcodes_df.empty:
            self.ambiguous_barcodes_df["proportion"] = pd.Series(dtype=float)
            self.ambiguous_barcodes_df["proportion_all_reads"] = pd.Series(
                dtype=float
            )
        else:
            # Rank non-exact/multi observations within each candidate profile.
            for profile, profile_df in self.ambiguous_barcodes_df.groupby("profile"):
                for prefix in ("Unique", "Multi"):
                    observed = profile_df["bc_name"] == prefix
                    sequence_counts = (
                        profile_df.loc[observed]
                        .groupby("bc_seq")["bc_count"]
                        .sum()
                        .to_dict()
                    )
                    ranked_names = {
                        sequence: f"{prefix}.{rank}"
                        for rank, sequence in enumerate(
                            sorted(
                                sequence_counts,
                                key=lambda sequence: (
                                    -sequence_counts[sequence],
                                    sequence,
                                ),
                            ),
                            start=1,
                        )
                    }
                    profile_observed = (
                        self.ambiguous_barcodes_df["profile"] == profile
                    ) & (self.ambiguous_barcodes_df["bc_name"] == prefix)
                    self.ambiguous_barcodes_df.loc[
                        profile_observed, "bc_name"
                    ] = self.ambiguous_barcodes_df.loc[
                        profile_observed, "bc_seq"
                    ].map(ranked_names)

            ambiguous_total = int(self.ambiguous_barcodes_df["bc_count"].sum())
            self.ambiguous_barcodes_df["proportion"] = (
                self.ambiguous_barcodes_df["bc_count"] / ambiguous_total
            )
            self.ambiguous_barcodes_df["proportion_all_reads"] = (
                0.0
                if total_reads == 0
                else self.ambiguous_barcodes_df["bc_count"] / total_reads
            )
        self.ambiguous_barcodes_df = self.ambiguous_barcodes_df.loc[
            :, ambiguous_columns
        ].sort_values(
            ["ambiguity_type", "bc_count", "bc_seq"],
            ascending=[True, False, True],
        ).reset_index(drop=True)

    def _build_output_frames(self) -> None:
        """Build summary, composition, and long barcode output tables."""
        total = self.outcome_counts["total_reads"]

        # Terminal outcomes are mutually exclusive and sum to total reads.
        terminal_outcomes = (
            "profile_not_found",
            "quality_filtered",
            "profile_ambiguous",
            "stock_unknown",
            "stock_ambiguous",
            "stock_assigned",
        )
        summary_rows = [("total_reads", total, 1.0 if total else 0.0)]
        summary_rows.extend(
            (
                outcome,
                self.outcome_counts[outcome],
                self._proportion(self.outcome_counts[outcome], total),
            )
            for outcome in terminal_outcomes
        )
        self.summary_df = pd.DataFrame(
            summary_rows,
            columns=["outcome", "read_count", "proportion_all_reads"],
        )

        # Profile table: denominator includes every uniquely classified profile.
        profile_total = sum(self.profile_counts.values())
        self.profiles_df = pd.DataFrame(
            [
                {
                    "profile": profile,
                    "index_structure": structure,
                    "read_count": count,
                    "proportion_all_reads": self._proportion(count, total),
                    "proportion_profile_classified": self._proportion(
                        count, profile_total
                    ),
                }
                for (profile, structure), count in sorted(self.profile_counts.items())
            ],
            columns=[
                "profile",
                "index_structure",
                "read_count",
                "proportion_all_reads",
                "proportion_profile_classified",
            ],
        )

        # Stock table: denominator includes every uniquely assigned stock/control.
        stock_total = sum(self.stock_counts.values())
        self.stocks_df = pd.DataFrame(
            [
                {
                    "profile": profile,
                    "stock": stock,
                    "ldist": ldist,
                    "read_count": count,
                    "proportion_all_reads": self._proportion(count, total),
                    "proportion_stock_assigned": self._proportion(count, stock_total),
                }
                for (profile, stock, ldist), count in sorted(
                    self.stock_counts.items()
                )
            ],
            columns=[
                "profile",
                "stock",
                "ldist",
                "read_count",
                "proportion_all_reads",
                "proportion_stock_assigned",
            ],
        )

        # Ambiguity table: retain profile ties and equal-nearest stock ties.
        self.ambiguous_df = pd.DataFrame(
            [
                {
                    "ambiguity_type": kind,
                    "labels": labels,
                    "ldist": ldist,
                    "read_count": count,
                    "proportion_all_reads": self._proportion(count, total),
                }
                for (kind, labels, ldist), count in sorted(
                    self.ambiguous_counts.items(),
                    key=lambda item: (
                        item[0][0],
                        item[0][1],
                        -1 if item[0][2] is None else item[0][2],
                    ),
                )
            ],
            columns=[
                "ambiguity_type",
                "labels",
                "ldist",
                "read_count",
                "proportion_all_reads",
            ],
        )

        self._build_barcode_frames(total_reads=total)

    @property
    def output_stem(self) -> str:
        """Return FASTQ basename with common FASTQ suffixes removed."""
        name = Path(self.settings.sample_path).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.lower().endswith(suffix):
                return name[: -len(suffix)]
        return Path(name).stem

    def settings_rows(self) -> list[tuple[str, object]]:
        """Return runtime settings plus every stock/control reference path."""
        rows: list[tuple[str, object]] = list(self.runtime_settings.items())

        # Preserve exact reference provenance in the emitted settings sheet.
        for profile, stocks in self.stock_reference_paths.items():
            for stock, path in stocks.items():
                rows.append((f"reference:{profile}:{stock}", path))
        return rows


# %% versions
"""
v1.3.0 20261110
    - Reports inexact reference matches with Parse-style edit suffixes.
    - Names single inexact matches Unique.N and tied matches Multi.N.

v1.2.0 20261110
    - Aggregates extracted barcodes before stock/reference matching.
    - Runs exact or distance matching once per distinct profile/barcode pair.
    - Reports barcode-reference matching progress with a determinate counter.

v1.1.0 20261110
    - Added parse-style long barcode tables for each classified profile.
    - Retains FASTA barcode names and slash-separated ambiguous stock labels.
    - Added barcode-level output for profile and stock ambiguities.

v1.0.0 20261110
    - Added standalone run-level FASTQ profile/stock screening.
    - Reuses parse-mode FASTQ streaming and structural extraction.
    - Supports exact or optional nearest-reference barcode assignment.
    - Builds workbook-ready summary, profile, stock, and ambiguity tables.
"""
