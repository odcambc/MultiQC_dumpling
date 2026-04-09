from pathlib import Path

import pandas as pd
import pytest

from multiqc import config
from multiqc_dumpling.dumpling_initialization import dumpling_plugin_execution_start
from multiqc_dumpling.modules.dumpling.dumpling import MultiqcModule


def make_module():
    return MultiqcModule.__new__(MultiqcModule)


@pytest.fixture(autouse=True)
def reset_multiqc_config():
    original_sp = dict(config.sp)
    original_ignore = list(getattr(config, "fn_ignore_files", []))
    original_plugin_config = getattr(config, "multiqc_dumpling", None)
    original_orf_length = getattr(config, "orf_length", None)
    had_n_variants = hasattr(config, "n_variants")
    original_n_variants = getattr(config, "n_variants", None)

    yield

    config.sp.clear()
    config.sp.update(original_sp)
    config.fn_ignore_files = original_ignore

    if original_plugin_config is None and hasattr(config, "multiqc_dumpling"):
        delattr(config, "multiqc_dumpling")
    else:
        config.multiqc_dumpling = original_plugin_config

    if original_orf_length is None and hasattr(config, "orf_length"):
        delattr(config, "orf_length")
    else:
        config.orf_length = original_orf_length

    if had_n_variants:
        config.n_variants = original_n_variants
    elif hasattr(config, "n_variants"):
        delattr(config, "n_variants")


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


def test_config_hook_requires_multiqc_dumpling_section():
    if hasattr(config, "multiqc_dumpling"):
        delattr(config, "multiqc_dumpling")

    with pytest.raises(ValueError, match="multiqc_dumpling"):
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
