"""Integration tests against real dumpling outputs.

Fixtures in ``tests/fixtures/baseline_run/`` are copies of files produced
by an actual dumpling pipeline run (see ``CLAUDE.md`` for the pipeline
this plugin consumes). They live in this repo so CI doesn't need a
checked-out dumpling tree — but their contract is "this is what dumpling
emits today." When dumpling changes a format, these tests should fail
loudly rather than silently producing empty plots.
"""

import fnmatch
import re
from pathlib import Path

import pytest

from multiqc import config
from multiqc_dumpling.dumpling_initialization import dumpling_plugin_execution_start
from multiqc_dumpling.modules.dumpling.dumpling import MultiqcModule


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "baseline_run"
SCORING_FIXTURE_DIR = Path(__file__).parent / "fixtures" / "scoring"


def make_module():
    return MultiqcModule.__new__(MultiqcModule)


def _read(filename: str) -> str:
    return (FIXTURE_DIR / filename).read_text()


# ──────────────────────────────────────────────────────────────
# Parser tests against real outputs
# ──────────────────────────────────────────────────────────────


def test_real_total_stats_parses():
    """Pipeline writes ``key\\tvalue`` lines (no header). Make sure we agree."""
    stats = MultiqcModule._parse_kv_stats_file(_read("baseline_1_total_processing.tsv"))
    assert stats == {
        "total_counts": 41122,
        "total_rejected_counts": 3589,
        "total_accepted_counts": 37533,
    }
    # Sanity: pipeline invariant.
    assert stats["total_accepted_counts"] + stats["total_rejected_counts"] == stats["total_counts"]


def test_real_accepted_stats_parses():
    stats = MultiqcModule._parse_kv_stats_file(_read("baseline_1_accepted_processing.tsv"))
    assert stats["accepted_syn_counts"] == 3726
    assert stats["accepted_sub_counts"] == 32242
    assert stats["accepted_stop_counts"] == 1565


def test_real_rejected_stats_parses():
    stats = MultiqcModule._parse_kv_stats_file(_read("baseline_1_rejected_processing.tsv"))
    assert stats["outside_orf_counts"] == 831
    assert stats["fs_counts"] == 545
    assert stats["wrong_variant_counts"] == 2213


def test_real_bbduk_decontam_log_parses():
    stats = MultiqcModule._parse_bbduk_log(_read("baseline_1.clean.bbduk.log"))
    assert stats["reads_in"] == 199_788
    assert stats["reads_out"] == 199_782


def test_real_bbduk_trim_log_parses():
    stats = MultiqcModule._parse_bbduk_log(_read("baseline_1.trim.bbduk.log"))
    assert stats["reads_in"] == 200_000
    # Trim removes a small number; just assert it actually dropped some.
    assert stats["reads_out"] < stats["reads_in"]


def test_real_bbmerge_log_parses():
    stats = MultiqcModule._parse_bbmerge_log(_read("baseline_1.bbmerge.log"))
    assert stats["pairs_in"] == 99_891
    assert stats["reads_in"] == 99_891 * 2
    assert stats["reads_out"] == 70_331
    assert stats["pct_joined"] == pytest.approx(70.408)


def test_real_bbmap_log_parses():
    """Real BBMap output uses ``Mapped reads:`` and a separate ``Percent mapped:``
    line, not the ``Mapped: N P%`` form. The parser needs to handle both."""
    stats = MultiqcModule._parse_bbmap_log(_read("baseline_1.bbmap_map.log"))
    assert stats["reads_in"] == 199_782
    assert stats["reads_out"] == 198_195
    assert stats["pct_mapped"] == pytest.approx(99.206, abs=0.01)


def test_scoring_backend_detection_from_path():
    backend = MultiqcModule._scoring_backend_from_path("/some/path/rosace")
    assert backend == "rosace"
    backend = MultiqcModule._scoring_backend_from_path("/some/path/lilace")
    assert backend == "lilace"
    backend = MultiqcModule._scoring_backend_from_path("/elsewhere/results")
    assert backend is None


def test_enrich2_replicate_name_from_root():
    n = MultiqcModule._enrich2_replicate_name_from_root
    assert n("/x/enrich/tsv/cond_A_R1_sel") == "cond_A_R1"
    assert n("/x/enrich/tsv/cond_B_R2_sel") == "cond_B_R2"
    assert n("/x/enrich/tsv/example_experiment_exp") is None  # _exp, not _sel
    assert n("/elsewhere") is None


def test_real_enrich2_scores_parses(monkeypatch):
    """Enrich2 per-replicate scores: 2372 variants, 103 significant at p<0.05.
    Pinned against the fixture; if enrich2 changes its column names, this fails."""
    monkeypatch.setattr(config, "n_variants", 2372, raising=False)

    stats, histogram, scores = make_module().parse_enrich2_scores(
        SCORING_FIXTURE_DIR / "enrich2_cond_A_R1_main_identifiers_scores.tsv"
    )

    assert stats["n_total"] == 2372
    assert stats["n_scored"] == 2372         # all rows have score + finite SE
    assert stats["pct_scored"] == pytest.approx(1.0)
    assert stats["n_significant"] == 103     # pvalue_raw < 0.05
    assert stats["fraction_significant"] == pytest.approx(103 / 2372, abs=0.001)
    assert len(histogram) > 0
    assert len(scores) == 2372


def test_parse_enrich2_missing_columns_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("hgvs\tscore\np.(A1V)\t0.5\n")  # missing SE, pvalue_raw
    with pytest.raises(ValueError, match="missing required columns"):
        make_module().parse_enrich2_scores(bad)


def test_real_rosace_scores_parses(monkeypatch):
    """Rosace's cond_A_scores.csv: 1770 variants, 400 significant at lfsr<0.05.
    Pinned against the fixture; if rosace changes its format, this fails loudly."""
    monkeypatch.setattr(config, "n_variants", 1770, raising=False)

    stats, histogram, scores = make_module().parse_scores(
        SCORING_FIXTURE_DIR / "cond_A_scores.csv"
    )

    assert stats["n_total"] == 1770
    assert stats["n_scored"] == 1770          # all rows in this file have mean+sd
    assert stats["pct_scored"] == pytest.approx(1.0)
    assert stats["n_significant"] == 400      # lfsr < 0.05 cutoff
    assert stats["fraction_significant"] == pytest.approx(400 / 1770, abs=0.001)
    assert stats["mean_score"] < 0            # selection screen → mostly negative means
    assert len(histogram) > 0
    assert len(scores) == 1770


def test_parse_scores_missing_columns_raises(tmp_path):
    """Format drift should fail loudly, not silently produce an empty result."""
    bad = tmp_path / "bad_scores.csv"
    bad.write_text("variants,position\np.(A1V),1\n")
    with pytest.raises(ValueError, match="missing required columns"):
        make_module().parse_scores(bad)


def test_parse_scores_excludes_inf_sd(tmp_path):
    """SD = +inf should not count as 'scored' per the user spec."""
    bad = tmp_path / "scores.csv"
    bad.write_text(
        "variants,position,mean,sd,lfsr\n"
        "p.(A1V),1,1.0,inf,0.01\n"
        "p.(A1L),1,0.5,0.1,0.4\n"
    )
    stats, _, _ = make_module().parse_scores(bad)
    assert stats["n_total"] == 2
    assert stats["n_scored"] == 1   # the inf-SD row is excluded


def test_real_samtools_stats_parses():
    """samtools stats is the only mapping source available in minimap2 runs."""
    stats = MultiqcModule._parse_samtools_stats(_read("baseline_1_samtools_stats.txt"))
    assert stats["reads_in"] == 353_908
    assert stats["reads_out"] == 350_389
    assert stats["pct_mapped"] == pytest.approx(99.005, abs=0.01)
    # Base counts for the base-survival waterfall.
    assert stats["bases_in"] == 106_160_357
    assert stats["bases_out"] == 103_245_630


def test_real_bbduk_log_extracts_bases():
    stats = MultiqcModule._parse_bbduk_log(_read("baseline_1.clean.bbduk.log"))
    assert stats["bases_in"] == 59_930_232
    assert stats["bases_out"] == 59_928_430


def test_real_bbmap_log_extracts_bases():
    stats = MultiqcModule._parse_bbmap_log(_read("baseline_1.bbmap_map.log"))
    # "Reads Used: N (M bases)" gives input bases; "Mapped bases: K" the output.
    assert "bases_in" in stats
    assert "bases_out" in stats
    assert stats["bases_in"] > stats["bases_out"] > 0


def test_extract_bbtools_version_from_real_log():
    """BBTools logs all start with ``Version N.NN`` after the executing line."""
    v = MultiqcModule._extract_tool_version(_read("baseline_1.clean.bbduk.log"), "BBTools")
    assert v == "39.13"


def test_extract_samtools_version_from_real_stats():
    v = MultiqcModule._extract_tool_version(_read("baseline_1_samtools_stats.txt"), "samtools")
    assert v == "1.23.1+htslib-1.23.1"


def test_extract_unknown_tool_returns_none():
    v = MultiqcModule._extract_tool_version("# Some random log content\n", "nonexistent_tool")
    assert v is None


def test_real_counts_file_parses(monkeypatch):
    """The counts CSV from process_counts.py should parse end-to-end."""
    monkeypatch.setattr(config, "n_variants", 4287, raising=False)
    plot_data, stats = make_module().parse_counts(FIXTURE_DIR / "baseline_1.csv")
    assert len(plot_data) > 0
    assert stats["Mean counts"] > 0
    assert stats["Median counts"] >= 0
    # The HGVS classifier should hit at least one missense from real data.
    assert stats.get("n_missense", 0) > 0


def test_real_refcoverage_parses(monkeypatch):
    monkeypatch.setattr(config, "n_variants", 4287, raising=False)
    plot_data, stats = make_module().parse_coverage(FIXTURE_DIR / "baseline_1.refCoverage")
    assert len(plot_data) > 0
    assert stats["Mean coverage"] > 0
    assert stats["Max coverage"] >= stats["Min coverage"]


# ──────────────────────────────────────────────────────────────
# Search-pattern coverage tests
# ──────────────────────────────────────────────────────────────


def _filename_matches_sp(filename: str, sp_entry: dict) -> bool:
    """Check if ``filename`` would match the search pattern entry.

    MultiQC supports either ``fn`` (fnmatch glob) or ``fn_re`` (regex) — both
    as strings. The test mirrors that logic so we verify what MultiQC itself
    will do at run time.
    """
    fn = sp_entry.get("fn")
    if isinstance(fn, str) and fnmatch.fnmatch(filename, fn):
        return True
    fn_re = sp_entry.get("fn_re")
    if isinstance(fn_re, str) and re.search(fn_re, filename):
        return True
    return False


@pytest.mark.parametrize(
    "sp_key, fixture_file",
    [
        ("dumpling/bbduk_decontam", "baseline_1.clean.bbduk.log"),
        ("dumpling/bbduk_trim",     "baseline_1.trim.bbduk.log"),
        ("dumpling/bbmerge",        "baseline_1.bbmerge.log"),
        ("dumpling/bbmap",          "baseline_1.bbmap_map.log"),
        ("dumpling/samtools_stats", "baseline_1_samtools_stats.txt"),
        ("dumpling/filter_stats",   "baseline_1_total_processing.tsv"),
        ("dumpling/accepted_stats", "baseline_1_accepted_processing.tsv"),
        ("dumpling/rejected_stats", "baseline_1_rejected_processing.tsv"),
        ("dumpling/counts",         "baseline_1.csv"),
    ],
)
def test_search_pattern_matches_real_filename(sp_key, fixture_file, _ensure_patterns_registered):
    """Each search pattern must match the filename dumpling actually emits.

    If this fails, the plugin's waterfall / counts / stats section will be silently empty
    for real pipeline runs.
    """
    sp_entry = config.sp.get(sp_key, {})
    assert sp_entry, f"No search-pattern entry registered for {sp_key}"
    assert _filename_matches_sp(fixture_file, sp_entry), (
        f"Search pattern {sp_key} ({sp_entry}) did not match real filename {fixture_file!r}."
    )


@pytest.fixture
def _ensure_patterns_registered():
    """Register the plugin's search patterns into config.sp for the test."""
    # The hook also wants to set ORF/n_variants etc., so seed a minimal config.
    variants_csv = FIXTURE_DIR.parent / "_tmp_variants.csv"
    variants_csv.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")
    config.multiqc_dumpling = {"orf": "1-9", "variants_file": str(variants_csv)}
    try:
        dumpling_plugin_execution_start()
        yield
    finally:
        variants_csv.unlink(missing_ok=True)
        if hasattr(config, "multiqc_dumpling"):
            delattr(config, "multiqc_dumpling")
