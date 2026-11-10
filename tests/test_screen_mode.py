from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from bcparse.config import (
    base_opts,
    dual_index_opts,
    get_screen_stock_settings,
    screen_stock_exclusions,
    settings_opts,
)
from bcparse.emit.workbooks import write_screen_workbook
from bcparse.modes.parse_mode import (
    fq_profile_parser,
    stream_fastq,
    stream_to_countdict,
)
from bcparse.modes.screen_mode import ScreenMode, StockMatch
from bcparse.settings import ScreenSettings


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(table)[::-1]


def parser_settings(profile: str, *, dual: bool = False, mismatches: int = 0) -> dict:
    settings = base_opts[profile].copy()
    settings.update(
        base_profile=profile,
        dualindex=False,
        mismatches=mismatches,
        mask=False,
        mask_qual=20,
    )
    if dual:
        settings.update(dual_index_opts[profile])
    return settings


def single_read(profile: str, barcode: str, *, p5_anchor: str | None = None) -> str:
    settings = base_opts[profile]
    ref_p5 = settings["ref_p5"] if p5_anchor is None else p5_anchor
    if settings["rdir"] == "fwd":
        return "ACGTACGT" + ref_p5 + barcode

    if settings["tdir_bc"] == "upstream":
        oriented_target = barcode + settings["ref_bc"]
    else:
        oriented_target = settings["ref_bc"] + barcode
    return "ACGTACGT" + ref_p5 + reverse_complement(oriented_target)


@pytest.mark.parametrize(
    ("profile", "barcode"),
    [
        ("M", "A" * 34),
        ("X", "C" * 45),
        ("SL8_bc", "G" * 34),
        ("SL8_epi", "T" * 24),
        ("pre_VIF", "AC" * 17),
        ("post_VIF", "GT" * 17),
    ],
)
def test_single_index_parser_extracts_every_screen_profile(
    profile: str,
    barcode: str,
) -> None:
    parser = fq_profile_parser(parser_settings(profile))
    sequence = single_read(profile, barcode)
    result = parser.parse(seq=sequence, qual="I" * len(sequence))

    assert result.bc_seq == barcode
    assert result.p5_seq == "ACGTACGT"


def test_parser_allows_fuzzy_structural_anchors_but_not_beyond_limit() -> None:
    profile = "X"
    barcode = "C" * 45
    ref_p5 = base_opts[profile]["ref_p5"]
    one_substitution = ("C" if ref_p5[0] != "C" else "A") + ref_p5[1:]
    sequence = single_read(profile, barcode, p5_anchor=one_substitution)

    fuzzy = fq_profile_parser(parser_settings(profile, mismatches=1))
    exact = fq_profile_parser(parser_settings(profile, mismatches=0))

    assert fuzzy.parse(seq=sequence, qual="I" * len(sequence)).bc_seq == barcode
    assert not hasattr(exact.parse(seq=sequence, qual="I" * len(sequence)), "bc_seq")


def test_dual_m_parser_extracts_and_repairs_short_barcode() -> None:
    settings = parser_settings("M", dual=True)
    parser = fq_profile_parser(settings)
    short_barcode = "G" * 30
    expected = settings["fill_seq"][-4:] + short_barcode
    oriented_target = (
        "TGCATGCA" + settings["ref_p7"] + short_barcode + settings["ref_bc"]
    )
    sequence = (
        "ACGTACGT" + settings["ref_p5"] + reverse_complement(oriented_target)
    )

    result = parser.parse(seq=sequence, qual="I" * len(sequence))

    assert result.p7_seq == "TGCATGCA"
    assert result.bc_seq == expected


def test_normal_parse_count_stream_uses_shared_profile_parser(tmp_path: Path) -> None:
    barcode = "A" * 34
    sequence = single_read("M", barcode)
    fastq_path = tmp_path / "parse.fastq"
    write_fastq(fastq_path, [("read", sequence, "I" * len(sequence))])
    settings = parser_settings("M")
    settings.update(sample_path=str(fastq_path), mean_qual=30)

    assert stream_to_countdict(settings) == {"ACGTACGT": {barcode: 1}}


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in records.items())
    )


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    path.write_text(
        "".join(
            f"@{name}\n{sequence}\n+\n{quality}\n"
            for name, sequence, quality in records
        )
    )


def test_stream_fastq_reads_gzipped_fastq(tmp_path: Path) -> None:
    fastq_path = tmp_path / "screen.fastq.gz"
    with gzip.open(fastq_path, "wt") as fastq:
        fastq.write("@read\nACGT\n+\nIIII\n")

    assert list(stream_fastq(str(fastq_path))) == [
        ("@read", "ACGT", "+", "IIII")
    ]


def test_screen_mode_terminal_outcomes_and_workbook(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    barcode_unique = "A" * 34
    barcode_ambiguous = "C" * 34
    barcode_unknown = "G" * 34
    barcode_near_unique = "C" + barcode_unique[1:]
    barcode_fuzzy_tie = barcode_unique[:-1] + "C"
    stock_b_tie = barcode_unique[:-2] + "CC"
    stock_a = tmp_path / "stock_a.fasta"
    stock_b = tmp_path / "stock_b.fasta"
    write_fasta(
        stock_a,
        {
            "A.1": barcode_unique,
            "A.2": barcode_ambiguous,
            "A.3": stock_b_tie,
        },
    )
    write_fasta(stock_b, {"B.1": barcode_ambiguous})

    monkeypatch.setattr(
        "bcparse.modes.screen_mode.get_screen_stock_settings",
        lambda: {
            "M": {"StockA": str(stock_a), "StockB": str(stock_b)},
            "X": {},
            "SL8_bc": {},
            "SL8_epi": {},
            "pre_VIF": {},
            "post_VIF": {},
        },
    )

    match_calls: dict[tuple[str, str], int] = {}
    original_match_stocks = ScreenMode._match_stocks

    def counted_match_stocks(
        self: ScreenMode,
        *,
        profile: str,
        barcode: str,
    ) -> StockMatch:
        key = (profile, barcode)
        match_calls[key] = match_calls.get(key, 0) + 1
        return original_match_stocks(self, profile=profile, barcode=barcode)

    monkeypatch.setattr(ScreenMode, "_match_stocks", counted_match_stocks)

    sl8_settings = base_opts["SL8_bc"]
    sl8_oriented = (
        "A" * 34
        + sl8_settings["ref_bc"]
        + base_opts["SL8_epi"]["ref_bc"]
        + "T" * 24
    )
    ambiguous_profile = (
        "ACGTACGT"
        + sl8_settings["ref_p5"]
        + reverse_complement(sl8_oriented)
    )

    sequences = {
        "assigned": single_read("M", barcode_unique),
        "assigned_repeat": single_read("M", barcode_unique),
        "stock_ambiguous": single_read("M", barcode_ambiguous),
        "stock_unknown": single_read("M", barcode_unknown),
        "stock_nearest": single_read("M", barcode_near_unique),
        "stock_nearest_tie": single_read("M", barcode_fuzzy_tie),
        "profile_not_found": "ACGT" * 30,
        "quality_filtered": single_read("M", barcode_unique),
        "profile_ambiguous": ambiguous_profile,
    }
    fastq_path = tmp_path / "run.fastq"
    write_fastq(
        fastq_path,
        [
            (
                name,
                sequence,
                "!" * len(sequence) if name == "quality_filtered" else "I" * len(sequence),
            )
            for name, sequence in sequences.items()
        ],
    )

    screen = ScreenMode(
        settings=ScreenSettings(
            sample_path=str(fastq_path),
            out_path=str(tmp_path),
            mismatches=0,
            dist_threshold=1,
        )
    )

    assert screen.outcome_counts["total_reads"] == 9
    assert screen.outcome_counts["stock_assigned"] == 4
    assert screen.outcome_counts["stock_ambiguous"] == 1
    assert screen.outcome_counts["stock_unknown"] == 1
    assert screen.outcome_counts["profile_not_found"] == 1
    assert screen.outcome_counts["quality_filtered"] == 1
    assert screen.outcome_counts["profile_ambiguous"] == 1
    assert sum(
        screen.outcome_counts[outcome]
        for outcome in (
            "profile_not_found",
            "quality_filtered",
            "profile_ambiguous",
            "stock_unknown",
            "stock_ambiguous",
            "stock_assigned",
        )
    ) == screen.outcome_counts["total_reads"]
    assert match_calls[("M", barcode_unique)] == 1
    assert set(match_calls.values()) == {1}
    assert set(screen.stocks_df["ldist"]) == {0, 1}
    assert set(
        screen.ambiguous_df.loc[
            screen.ambiguous_df["ambiguity_type"] == "stock",
            "ldist",
        ]
    ) == {0}

    m_barcodes = screen.barcode_dfs["M"]
    assert m_barcodes["bc_count"].sum() == 6
    assert m_barcodes.loc[
        m_barcodes["bc_seq"] == barcode_unique,
        "bc_count",
    ].iloc[0] == 2
    assert m_barcodes["proportion"].sum() == pytest.approx(1.0)
    assert "A.1" in set(m_barcodes["bc_name"])
    assert m_barcodes["bc_name"].str.startswith("Unique.").any()
    assert m_barcodes["bc_name"].str.startswith("Multi.").any()
    assert set(m_barcodes["stock"]) >= {
        "StockA",
        "StockA/StockB",
        "Unknown",
    }
    assert m_barcodes.loc[
        m_barcodes["bc_seq"] == barcode_near_unique,
        "ref_match",
    ].iloc[0] == "A.1_1C"
    assert m_barcodes.loc[
        m_barcodes["bc_seq"] == barcode_fuzzy_tie,
        "ref_match",
    ].iloc[0] == "A.1_34C/A.3_33A"
    assert barcode_fuzzy_tie not in set(screen.ambiguous_barcodes_df["bc_seq"])
    assert set(
        screen.ambiguous_barcodes_df.loc[
            screen.ambiguous_barcodes_df["ambiguity_type"] == "stock",
            "stocks",
        ]
    ) == {"StockA/StockB"}
    assert set(
        screen.ambiguous_barcodes_df.loc[
            screen.ambiguous_barcodes_df["ambiguity_type"] == "profile",
            "profiles",
        ]
    ) == {"SL8_bc/SL8_epi"}

    workbook_path = write_screen_workbook(screen)
    assert Path(workbook_path).name == "run_Screen.xlsx"
    assert pd.ExcelFile(workbook_path).sheet_names == [
        "Summary",
        "M barcodes",
        "Ambiguous barcodes",
        "Screen settings",
    ]
    summary_sheet = pd.read_excel(workbook_path, sheet_name="Summary", header=None)
    assert {"Read outcomes", "Profiles", "Stocks", "Ambiguous"}.issubset(
        set(summary_sheet[0].dropna())
    )


def test_screen_stock_settings_follow_config_with_only_explicit_exclusions() -> None:
    actual = get_screen_stock_settings()
    expected = {
        (stock_settings["base_profile"], stock)
        for stock, stock_settings in settings_opts.items()
        if stock not in screen_stock_exclusions
        and stock_settings.get("barcode_path")
    }
    expected.update((profile, "Spike") for profile in base_opts)

    assert {
        (profile, stock)
        for profile, stocks_by_profile in actual.items()
        for stock in stocks_by_profile
    } == expected
    assert not (screen_stock_exclusions & actual["M"].keys())


def test_screen_dist_threshold_defaults_to_exact_matching() -> None:
    assert ScreenSettings().dist_threshold == 0
