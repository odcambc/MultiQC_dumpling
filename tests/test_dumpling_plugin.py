from pathlib import Path

import pandas as pd
import pytest

from multiqc import config
from multiqc_dumpling.dumpling_initialization import dumpling_plugin_execution_start
from multiqc_dumpling.modules.dumpling.dumpling import MultiqcModule


def make_module():
    return MultiqcModule.__new__(MultiqcModule)


_HOOK_SET_ATTRS = [
    "orf_length", "orf_start", "orf_end", "orf_n_codons",
    "n_variants", "dumpling_settings",
    "amplicon_length", "amplicon_start", "amplicon_end", "minq", "min_variants",
    "dumpling_pipeline_version", "dumpling_pipeline_config",
    "sample_names_replace",
]


@pytest.fixture(autouse=True)
def reset_multiqc_config():
    original_sp = dict(config.sp)
    original_ignore = list(getattr(config, "fn_ignore_files", []))
    original_clean_exts = list(getattr(config, "fn_clean_exts", []) or [])
    original_plugin_config = getattr(config, "multiqc_dumpling", None)
    saved = {attr: (hasattr(config, attr), getattr(config, attr, None)) for attr in _HOOK_SET_ATTRS}

    yield

    config.sp.clear()
    config.sp.update(original_sp)
    config.fn_ignore_files = original_ignore
    config.fn_clean_exts = original_clean_exts

    if original_plugin_config is None and hasattr(config, "multiqc_dumpling"):
        delattr(config, "multiqc_dumpling")
    else:
        config.multiqc_dumpling = original_plugin_config

    for attr, (had, val) in saved.items():
        if had:
            setattr(config, attr, val)
        elif hasattr(config, attr):
            delattr(config, attr)


def test_parse_counts_zero_counts_returns_plot_data_first(tmp_path):
    counts_path = tmp_path / "counts.csv"
    pd.DataFrame(
        {
            "mutation": ["a", "b"],
            "length": [1, 1],
            "hgvs": ["p.A1A", "p.A2A"],
            "count": [0, 0],
        }
    ).to_csv(counts_path, index=False)

    config.n_variants = 2
    plot_data, stats_data = make_module().parse_counts(counts_path)

    assert plot_data == {0: 2}
    assert stats_data["Number of zero counts"] == 2
    assert stats_data["Fraction of zero counts"] == 1


def test_parse_coverage_respects_requested_max_bins(tmp_path):
    coverage_path = tmp_path / "sample.refCoverage"
    pd.DataFrame({"RefPos": range(1, 101), "Coverage": range(1000, 1100)}).to_csv(
        coverage_path, sep="\t", index=False
    )

    config.n_variants = 1
    plot_data, stats_data = make_module().parse_coverage(coverage_path, bin_n=10)

    assert len(plot_data) == 10
    assert stats_data["Max coverage"] == 1099
    assert stats_data["Min coverage"] == 1000


def test_config_hook_returns_silently_without_multiqc_dumpling_section():
    if hasattr(config, "multiqc_dumpling"):
        delattr(config, "multiqc_dumpling")

    # Should not raise — absent section means this isn't a dumpling run.
    dumpling_plugin_execution_start()
    assert not hasattr(config, "orf_length")


def test_config_hook_raises_when_orf_missing():
    config.multiqc_dumpling = {"variants_file": "/nonexistent"}

    with pytest.raises(ValueError, match="orf"):
        dumpling_plugin_execution_start()


def test_config_hook_requires_variants_file_to_exist(tmp_path):
    config.multiqc_dumpling = {
        "orf": "1-10",
        "variants_file": str(tmp_path / "missing.csv"),
    }

    with pytest.raises(FileNotFoundError, match="variants_file"):
        dumpling_plugin_execution_start()


def test_config_hook_sets_n_variants_and_ignore_file(tmp_path):
    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n2,p.A2A\n", encoding="utf-8")

    config.multiqc_dumpling = {
        "orf": "5-13",
        "variants_file": str(variants_path),
    }

    dumpling_plugin_execution_start()

    assert config.orf_length == 9
    assert config.n_variants == 2
    assert str(variants_path) in config.fn_ignore_files


def test_config_hook_stores_orf_start_end_and_codons(tmp_path):
    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")

    config.multiqc_dumpling = {
        "orf": "4-12",
        "variants_file": str(variants_path),
    }

    dumpling_plugin_execution_start()

    assert config.orf_start == 4
    assert config.orf_end == 12
    assert config.orf_length == 9
    assert config.orf_n_codons == 3


def test_config_hook_reads_pipeline_version_and_config(tmp_path):
    """Pipeline reproducibility metadata is plumbed through cleanly when present."""
    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")

    pipeline_cfg = tmp_path / "dumpling_config.yaml"
    pipeline_cfg.write_text(
        "experiment: test_exp\norf: '141-1568'\nscoring_backend: rosace\n",
        encoding="utf-8",
    )

    config.multiqc_dumpling = {
        "orf": "1-10",
        "variants_file": str(variants_path),
        "pipeline_version": "v0.2.0-5-gdeadbee",
        "pipeline_config_file": str(pipeline_cfg),
    }

    dumpling_plugin_execution_start()

    assert config.dumpling_pipeline_version == "v0.2.0-5-gdeadbee"
    assert config.dumpling_pipeline_config["experiment"] == "test_exp"
    assert config.dumpling_pipeline_config["scoring_backend"] == "rosace"


def test_experiment_file_registers_unambiguous_renames(tmp_path):
    """Each fastq prefix used by a single sample should produce a rename rule."""
    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")

    exp_path = tmp_path / "exp.csv"
    exp_path.write_text(
        "sample,condition,replicate,time,tile,file\n"
        "sampA,c1,1,0,1,fastq_001\n"
        "sampB,c2,1,0,1,fastq_002\n",
        encoding="utf-8",
    )

    config.multiqc_dumpling = {
        "orf": "1-10",
        "variants_file": str(variants_path),
        "experiment_file": str(exp_path),
    }
    dumpling_plugin_execution_start()

    assert config.sample_names_replace.get("fastq_001") == "sampA"
    assert config.sample_names_replace.get("fastq_002") == "sampB"


def test_experiment_file_skips_ambiguous_prefixes(tmp_path, caplog):
    """When multiple samples share a fastq, we must NOT pick one arbitrarily."""
    import logging

    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")

    exp_path = tmp_path / "exp.csv"
    exp_path.write_text(
        "sample,condition,replicate,time,tile,file\n"
        "sampA,c1,1,0,1,shared_fastq\n"
        "sampB,c2,1,0,1,shared_fastq\n"
        "sampC,c1,1,0,1,unique_fastq\n",
        encoding="utf-8",
    )

    config.multiqc_dumpling = {
        "orf": "1-10",
        "variants_file": str(variants_path),
        "experiment_file": str(exp_path),
    }
    with caplog.at_level(logging.WARNING, logger="multiqc"):
        dumpling_plugin_execution_start()

    renames = getattr(config, "sample_names_replace", {}) or {}
    assert "shared_fastq" not in renames                  # ambiguous → skipped
    assert renames.get("unique_fastq") == "sampC"         # unambiguous → renamed
    assert any("shared by multiple samples" in r.message for r in caplog.records)


def test_config_hook_warns_on_missing_pipeline_config(tmp_path, caplog):
    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")

    config.multiqc_dumpling = {
        "orf": "1-10",
        "variants_file": str(variants_path),
        "pipeline_config_file": str(tmp_path / "missing.yaml"),
    }
    dumpling_plugin_execution_start()
    # Should not raise — but should log a warning, and pipeline_config stays None.
    assert config.dumpling_pipeline_config is None


def test_config_hook_reads_optional_settings(tmp_path):
    variants_path = tmp_path / "variants.csv"
    variants_path.write_text("pos,hgvs\n1,p.A1A\n", encoding="utf-8")

    config.multiqc_dumpling = {
        "orf": "1-10",
        "variants_file": str(variants_path),
        "amplicon_length": "500",
        "minq": "30",
        "min_variants": "2",
    }

    dumpling_plugin_execution_start()

    assert config.amplicon_length == 500
    assert config.minq == 30
    assert config.min_variants == 2
    assert config.dumpling_settings == {"amplicon_length": 500, "minq": 30, "min_variants": 2}


def test_classify_hgvs_missense():
    m = MultiqcModule._classify_hgvs
    assert m("p.A1V") == "missense"
    assert m("p.G10D") == "missense"


def test_classify_hgvs_synonymous():
    m = MultiqcModule._classify_hgvs
    assert m("p.A1A") == "synonymous"
    assert m("p.G10G") == "synonymous"


def test_classify_hgvs_nonsense():
    m = MultiqcModule._classify_hgvs
    assert m("p.A1*") == "nonsense"
    assert m("p.G10*") == "nonsense"


def test_classify_hgvs_indel_types():
    m = MultiqcModule._classify_hgvs
    assert m("p.A1del") == "deletion"
    assert m("p.A1insGGG") == "insertion"
    assert m("p.A1delinsGGG") == "indel"


def test_classify_hgvs_wildtype():
    m = MultiqcModule._classify_hgvs
    assert m("p.=") == "wildtype"
    assert m("WT") == "wildtype"


# ──────────────────────────────────────────────────────────────
# Pipeline step log parsers
# ──────────────────────────────────────────────────────────────

_BBDUK_DECONTAM_LOG = """\
Input:                  1000000 reads           150000000 bases.
Contaminants:           50000 reads (5.0000%)   7500000 bases (5.0000%)
Result:                 950000 reads (95.0000%)  142500000 bases (95.0000%)
"""

_BBDUK_TRIM_LOG = """\
Input:                  950000 reads            142500000 bases.
QTrimmed:               0 reads (0.0000%)       10000 bases (0.0007%)
KTrimmed:               20000 reads (2.1053%)   3000000 bases (2.1053%)
Total Removed:          20000 reads (2.1053%)   3010000 bases (2.1123%)
Result:                 930000 reads (97.8947%)  139490000 bases (97.8877%)
"""

_BBMERGE_LOG = """\
Pairs:                  465000
Joined:                 400000          86.02%
Ambiguous:              0                0.00%
No Solution:            65000           13.98%
Avg Insert:             300
"""

_BBMAP_LOG = """\
Reads Used:             400000          (120000000 bases)
Mapped:                 380000  95.0000%        (114000000 bases)
Unambiguous:            375000  93.7500%        (112500000 bases)
Ambiguous:              5000     1.2500%        (1500000 bases)
Low Quality:            0        0.0000%        (0 bases)
"""

# Dumpling's process_variants.write_stats_file emits long-format key\tvalue
# lines (no header). The same shape is used for total / accepted / rejected.
_TOTAL_STATS_TSV = (
    "total_counts\t380000\n"
    "total_rejected_counts\t80000\n"
    "total_accepted_counts\t300000\n"
)
_ACCEPTED_STATS_TSV = (
    "accepted_syn_counts\t3612\n"
    "accepted_sub_counts\t31491\n"
    "accepted_stop_counts\t1871\n"
    "accepted_ins_counts\t0\n"
    "accepted_del_counts\t0\n"
    "accepted_insdel_counts\t0\n"
)
_REJECTED_STATS_TSV = (
    "outside_orf_counts\t821\n"
    "fs_counts\t358\n"
    "wrong_codon_counts\t0\n"
    "wrong_variant_counts\t1865\n"
    "insdel_variant_counts\t0\n"
    "multi_variant_counts\t0\n"
    "unexpected_mutation_counts\t0\n"
)


def test_parse_bbduk_log_decontam():
    stats = MultiqcModule._parse_bbduk_log(_BBDUK_DECONTAM_LOG)
    assert stats["reads_in"] == 1_000_000
    assert stats["reads_out"] == 950_000
    assert stats["bases_in"] == 150_000_000


def test_parse_bbduk_log_trim():
    stats = MultiqcModule._parse_bbduk_log(_BBDUK_TRIM_LOG)
    assert stats["reads_in"] == 950_000
    assert stats["reads_out"] == 930_000


def test_parse_bbmerge_log():
    stats = MultiqcModule._parse_bbmerge_log(_BBMERGE_LOG)
    assert stats["pairs_in"] == 465_000
    assert stats["reads_in"] == 465_000 * 2   # pairs × 2
    assert stats["reads_out"] == 400_000
    assert stats["pct_joined"] == pytest.approx(86.02)


def test_parse_bbmerge_log_reads_format():
    """BBMerge alternate output style: 'Joined: N reads (P%)' instead of 'Joined: N  P%'."""
    log_content = (
        "Input:                  930000 reads           139490000 bases.\n"
        "Joined:                 800000 reads (86.02%)           240000000 bases.\n"
        "No Solution:            130000 reads (13.98%)           38900000 bases.\n"
        "Avg Insert:             300\n"
    )
    stats = MultiqcModule._parse_bbmerge_log(log_content)
    assert stats["reads_out"] == 800_000
    assert stats["pct_joined"] == pytest.approx(86.02)


def test_parse_bbmap_log():
    stats = MultiqcModule._parse_bbmap_log(_BBMAP_LOG)
    assert stats["reads_in"] == 400_000
    assert stats["reads_out"] == 380_000
    assert stats["pct_mapped"] == pytest.approx(95.0)


def test_parse_total_stats():
    stats = MultiqcModule._parse_kv_stats_file(_TOTAL_STATS_TSV)
    assert stats == {
        "total_counts": 380_000,
        "total_rejected_counts": 80_000,
        "total_accepted_counts": 300_000,
    }


def test_parse_accepted_stats():
    stats = MultiqcModule._parse_kv_stats_file(_ACCEPTED_STATS_TSV)
    assert stats["accepted_syn_counts"] == 3612
    assert stats["accepted_sub_counts"] == 31491
    assert stats["accepted_stop_counts"] == 1871


def test_parse_rejected_stats():
    stats = MultiqcModule._parse_kv_stats_file(_REJECTED_STATS_TSV)
    assert stats["outside_orf_counts"] == 821
    assert stats["wrong_variant_counts"] == 1865
    assert stats["fs_counts"] == 358


def test_parse_kv_stats_skips_blank_and_malformed_lines():
    content = "good_key\t42\n\nbad_line_no_tab\nanother_good\t7\n"
    stats = MultiqcModule._parse_kv_stats_file(content)
    assert stats == {"good_key": 42, "another_good": 7}


def test_compute_processing_metrics_derives_ratios():
    mod = make_module()
    mod.pipeline_filter_data = {
        "sample1": {"total_counts": 1000, "total_accepted_counts": 800, "total_rejected_counts": 200},
    }
    mod.pipeline_rejected_data = {
        "sample1": {"outside_orf_counts": 100, "fs_counts": 50, "wrong_variant_counts": 50},
    }
    mod.pipeline_accepted_data = {}

    metrics = mod._compute_processing_metrics()

    assert metrics["sample1"]["read_acceptance_rate"] == pytest.approx(0.8)
    assert metrics["sample1"]["pct_reads_outside_orf"] == pytest.approx(0.1)
    assert metrics["sample1"]["pct_reads_frameshift"] == pytest.approx(0.05)
    assert metrics["sample1"]["pct_reads_wrong_variant"] == pytest.approx(0.05)


def test_compute_processing_metrics_skips_sample_with_no_total():
    mod = make_module()
    mod.pipeline_filter_data = {}
    mod.pipeline_rejected_data = {"orphan": {"outside_orf_counts": 1}}
    mod.pipeline_accepted_data = {}

    assert mod._compute_processing_metrics() == {}


def test_parse_counts_includes_variant_type_fractions(tmp_path):
    counts_path = tmp_path / "counts.csv"
    pd.DataFrame(
        {
            "mutation": ["a", "b", "c", "d"],
            "length": [1, 1, 1, 1],
            "hgvs": ["p.A1V", "p.A1A", "p.A1*", "p.A1del"],
            "count": [10, 20, 5, 8],
        }
    ).to_csv(counts_path, index=False)

    config.n_variants = 4
    _, stats = make_module().parse_counts(counts_path)

    assert stats["n_missense"] == 1
    assert stats["n_synonymous"] == 1
    assert stats["n_nonsense"] == 1
    assert stats["n_deletion"] == 1
    assert stats["fraction_missense"] == pytest.approx(0.25)
    assert stats["fraction_synonymous"] == pytest.approx(0.25)
