#!/usr/bin/env python

from __future__ import print_function

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, bargraph
import importlib_metadata


# Initialise the main MultiQC logger
log = logging.getLogger(__name__)
log.info("Loaded MultiQC_dumpling v%s", importlib_metadata.version("multiqc_dumpling"))

# TODO: validate the report on realistic data. The current ``example_experiment``
# fixture is contrived (multiple samples share the same input FASTQ, adapter
# content is <0.5%, base composition plots look unusually flat). Several plots
# that depend on a typical Illumina read distribution may render funny on this
# data — verify them on a real DMS library with proper adapter content and
# per-sample FASTQs before treating any flat/odd-looking plot as a bug in this
# plugin vs. a data-quality artifact.
_VARIANT_TYPES = ["wildtype", "synonymous", "missense", "nonsense", "insertion", "deletion", "indel", "other"]

# Significance threshold for rosace/lilace lfsr — conventional 0.05.
# A variant is "significant" if lfsr < this value.
LFSR_THRESHOLD = 0.05
# Enrich2 doesn't produce an lfsr; the closest comparable column is
# ``pvalue_raw`` (Wald-style test of score ≠ 0). Same conventional cutoff.
ENRICH2_PVALUE_THRESHOLD = 0.05


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="Dumpling",
            target="dumpling",
            anchor="dumpling",
            href="https://github.com/odcambc/dumpling",
            info=" is a module for conducting quality checks of DMS libraries.",
        )

        # Find load the counts files.
        self.dumpling_count_data = dict()
        self.dumpling_count_plot_data = dict()

        for f in self.find_log_files("dumpling/counts", filehandles=True):
            parsed_counts = self.parse_counts(f["f"])
            self.dumpling_count_plot_data[f["s_name"]] = parsed_counts[0]
            self.dumpling_count_data[f["s_name"]] = parsed_counts[1]

        # Find load the coverage results from gatk asm results.
        self.dumpling_coverage_data = dict()
        self.dumpling_coverage_plot_data = dict()

        for f in self.find_log_files(
            "gatk/analyze_saturation_mutagenesis/refcoverage", filehandles=True
        ):
            (
                self.dumpling_coverage_plot_data[f["s_name"]],
                self.dumpling_coverage_data[f["s_name"]],
            ) = self.parse_coverage(f["f"])

        # Parse pipeline step logs for the read-loss waterfall
        self.pipeline_decontam_data = {}
        self.pipeline_trim_data = {}
        self.pipeline_merge_data = {}
        self.pipeline_map_data = {}
        self.pipeline_filter_data = {}
        self.pipeline_accepted_data = {}
        self.pipeline_rejected_data = {}

        # Tool versions detected from log files — first hit wins per tool.
        self.tool_versions: dict[str, str] = {}

        def _scan_versions(content: str, *tools: str) -> None:
            for tool in tools:
                if tool in self.tool_versions:
                    continue
                v = self._extract_tool_version(content, tool)
                if v:
                    self.tool_versions[tool] = v

        for f in self.find_log_files("dumpling/bbduk_decontam"):
            parsed = self._parse_bbduk_log(f["f"])
            if parsed:
                self.pipeline_decontam_data[f["s_name"]] = parsed
            _scan_versions(f["f"], "BBTools")

        for f in self.find_log_files("dumpling/bbduk_trim"):
            parsed = self._parse_bbduk_log(f["f"])
            if parsed:
                self.pipeline_trim_data[f["s_name"]] = parsed
            _scan_versions(f["f"], "BBTools")

        for f in self.find_log_files("dumpling/bbmerge"):
            parsed = self._parse_bbmerge_log(f["f"])
            if parsed:
                self.pipeline_merge_data[f["s_name"]] = parsed
            _scan_versions(f["f"], "BBTools")

        for f in self.find_log_files("dumpling/bbmap"):
            parsed = self._parse_bbmap_log(f["f"])
            if parsed:
                self.pipeline_map_data[f["s_name"]] = parsed
            _scan_versions(f["f"], "BBTools")

        # samtools stats wins over bbmap log when both exist — it's the
        # only mapping source available in minimap2 runs and gives the same
        # answer as bbmap on bbmap runs (both reflect the post-mapping BAM).
        for f in self.find_log_files("dumpling/samtools_stats"):
            parsed = self._parse_samtools_stats(f["f"])
            if parsed:
                self.pipeline_map_data[f["s_name"]] = parsed
            _scan_versions(f["f"], "samtools")

        for f in self.find_log_files("dumpling/filter_stats"):
            parsed = self._parse_kv_stats_file(f["f"])
            if parsed:
                self.pipeline_filter_data[f["s_name"]] = parsed

        for f in self.find_log_files("dumpling/accepted_stats"):
            parsed = self._parse_kv_stats_file(f["f"])
            if parsed:
                self.pipeline_accepted_data[f["s_name"]] = parsed

        for f in self.find_log_files("dumpling/rejected_stats"):
            parsed = self._parse_kv_stats_file(f["f"])
            if parsed:
                self.pipeline_rejected_data[f["s_name"]] = parsed

        # Scoring outputs (rosace / lilace) — same per-condition CSV format
        # for both backends; the backend is inferred from the parent dir name.
        # Each key is "backend:condition" so reports can group cleanly.
        self.scoring_stats = {}           # for general stats / table
        self.scoring_plot_data = {}       # for histograms
        self.scoring_score_vectors = {}   # raw mean series for concordance
        # filecontents=False — we read the CSV with pandas directly from disk
        # rather than slurping it through MultiQC's loader.
        for f in self.find_log_files("dumpling/scores", filecontents=False):
            backend = self._scoring_backend_from_path(f["root"])
            if backend is None:
                continue
            file_path = Path(f["root"]) / f["fn"]
            try:
                stats, hist, scores = self.parse_scores(file_path)
            except Exception as exc:
                log.warning("Failed to parse score file %s: %s", file_path, exc)
                continue
            key = f"{backend}:{f['s_name']}"
            self.scoring_stats[key] = stats
            self.scoring_plot_data[key] = hist
            self.scoring_score_vectors[key] = scores

        # Enrich2 per-replicate scores. The MultiQC-derived s_name is just the
        # generic filename ``main_identifiers_scores``; the per-replicate
        # identity comes from the parent directory (``cond_A_R1_sel`` etc.).
        for f in self.find_log_files("dumpling/enrich2_scores", filecontents=False):
            replicate = self._enrich2_replicate_name_from_root(f["root"])
            if replicate is None:
                continue  # not a per-replicate _sel directory
            file_path = Path(f["root"]) / f["fn"]
            try:
                stats, hist, scores = self.parse_enrich2_scores(file_path)
            except Exception as exc:
                log.warning("Failed to parse enrich2 scores %s: %s", file_path, exc)
                continue
            key = f"enrich2:{replicate}"
            self.scoring_stats[key] = stats
            self.scoring_plot_data[key] = hist
            self.scoring_score_vectors[key] = scores

        # Filter out samples matching ignored sample names
        self.dumpling_count_plot_data = self.ignore_samples(
            self.dumpling_count_plot_data
        )
        self.dumpling_count_data = self.ignore_samples(self.dumpling_count_data)

        self.dumpling_coverage_data = self.ignore_samples(self.dumpling_coverage_data)
        self.dumpling_coverage_plot_data = self.ignore_samples(
            self.dumpling_coverage_plot_data
        )
        self.pipeline_decontam_data = self.ignore_samples(self.pipeline_decontam_data)
        self.pipeline_trim_data = self.ignore_samples(self.pipeline_trim_data)
        self.pipeline_merge_data = self.ignore_samples(self.pipeline_merge_data)
        self.pipeline_map_data = self.ignore_samples(self.pipeline_map_data)
        self.pipeline_filter_data = self.ignore_samples(self.pipeline_filter_data)
        self.pipeline_accepted_data = self.ignore_samples(self.pipeline_accepted_data)
        self.pipeline_rejected_data = self.ignore_samples(self.pipeline_rejected_data)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.dumpling_count_plot_data) == 0:
            log.debug("Could not find any reports in %s", config.analysis_dir)
            raise ModuleNoSamplesFound

        log.info("Found %s processed count files", len(self.dumpling_count_plot_data))

        if len(self.dumpling_coverage_data) == 0:
            log.warning("Could not find any coverage reports in %s", config.analysis_dir)
        else:
            log.info("Found %s coverage report files", len(self.dumpling_coverage_data))

        # Write all the parsed report data to files
        self.write_data_file(self.dumpling_count_data, "multiqc_dumpling_counts")
        if self.dumpling_coverage_data:
            self.write_data_file(self.dumpling_coverage_data, "multiqc_dumpling_coverage")
        self.write_data_file(self.dumpling_count_plot_data, "multiqc_dumpling_counts_plot")
        if self.dumpling_coverage_plot_data:
            self.write_data_file(self.dumpling_coverage_plot_data, "multiqc_dumpling_coverage_plot")

        # Add a settings/ORF info section at the top
        self._add_settings_section()

        # Pipeline survival waterfall — combined reads + bases with a built-in
        # dataset switcher. Skipped automatically if no step logs were found.
        self._add_pipeline_waterfall_section()

        # Scoring sections (rosace / lilace)
        if self.scoring_stats:
            self.scoring_stats = self.ignore_samples(self.scoring_stats)
            self._add_scoring_general_stats()
            self._add_scoring_histogram_section()
            self._add_scoring_concordance_section()
            self.write_data_file(self.scoring_stats, "multiqc_dumpling_scoring")

        # Add to the general stats table
        counts_headers = {
            "Mean counts": {
                "title": "Mean variant counts",
                "description": "The mean number of counts for a variant in the sample.",
                "min": 0,
                "scale": "RdYlGn-rev",
            },
            "Median counts": {
                "title": "Median sequencing depth",
                "description": "The median number of counts for a variant in the sample.",
                "min": 0,
                "scale": "RdYlGn-rev",
            },
            "Number of zero counts": {
                "title": "Number of zero counts",
                "description": "The number of variants with zero counts in the sample (i.e., missing).",
                "min": 0,
                "scale": "RdYlGn-rev",
            },
            "Fraction of zero counts": {
                "title": "Fraction of zero counts",
                "description": "The fraction of variants with zero counts.",
                "min": 0,
                "scale": "RdYlGn-rev",
            },
            "fraction_missense": {
                "title": "Missense frac.",
                "description": "Fraction of designed variants classified as missense substitutions (by HGVS notation).",
                "min": 0,
                "max": 1,
                "scale": "Blues",
                "format": "{:.2f}",
            },
            "fraction_synonymous": {
                "title": "Synonymous frac.",
                "description": "Fraction of designed variants classified as synonymous substitutions (by HGVS notation).",
                "min": 0,
                "max": 1,
                "scale": "Greens",
                "format": "{:.2f}",
            },
            "fraction_nonsense": {
                "title": "Nonsense frac.",
                "description": "Fraction of designed variants classified as nonsense (stop codon) mutations (by HGVS notation).",
                "min": 0,
                "max": 1,
                "scale": "Reds",
                "format": "{:.2f}",
            },
        }

        self.general_stats_addcols(self.dumpling_count_data, counts_headers)

        coverage_headers = {
            "Max coverage": {
                "title": "Max sequencing depth",
                "description": "The maximum sequencing depth of a position in the reference.",
                "min": 0,
                "scale": "RdYlGn-rev",
            },
            "Min coverage": {
                "title": "Min sequencing depth",
                "description": "The minimum sequencing depth of a position in the reference.",
                "min": 0,
                "scale": "RdYlGn-rev",
            },
            "Mean coverage": {
                "title": "Mean sequencing depth",
                "description": "The mean sequencing depth of a position in the reference.",
                "min": 0,
                "scale": "RdYlGn-rev",
            }
        }

        if self.dumpling_coverage_data:
            self.general_stats_addcols(self.dumpling_coverage_data, coverage_headers)

        # Derived processing metrics — combine total/accepted/rejected counts
        # into per-sample ratios that surface library-quality issues at a glance
        # (off-target reads, frameshifts, off-design variants, etc.).
        processing_metrics = self._compute_processing_metrics()
        if processing_metrics:
            processing_headers = {
                "read_acceptance_rate": {
                    "title": "Read acceptance",
                    "description": (
                        "Fraction of mapped reads kept after variant filtering "
                        "(total_accepted_counts / total_counts)."
                    ),
                    "min": 0,
                    "max": 1,
                    "scale": "RdYlGn",
                    "format": "{:.2%}",
                },
                "pct_reads_outside_orf": {
                    "title": "% outside ORF",
                    "description": "Fraction of reads with no in-ORF variant; flags amplicon-design issues.",
                    "min": 0,
                    "max": 1,
                    "scale": "Reds",
                    "format": "{:.2%}",
                    "hidden": True,
                },
                "pct_reads_frameshift": {
                    "title": "% frameshift",
                    "description": "Fraction of reads with a frameshift; flags library/synthesis errors.",
                    "min": 0,
                    "max": 1,
                    "scale": "Reds",
                    "format": "{:.2%}",
                    "hidden": True,
                },
                "pct_reads_wrong_variant": {
                    "title": "% off-design",
                    "description": "Fraction of reads carrying a variant not in the designed library.",
                    "min": 0,
                    "max": 1,
                    "scale": "Reds",
                    "format": "{:.2%}",
                    "hidden": True,
                },
            }
            self.general_stats_addcols(processing_metrics, processing_headers)
            self.write_data_file(processing_metrics, "multiqc_dumpling_processing")

        # Plot the counts histogram

        # Create line plot for counts histogram
        pconfig = {
            "id": "dumpling_counts_hist",
            "title": "Variant count histogram",
            "ylab": "Number of variants",
            "xlab": "Counts of variant in sample",
            "ymin": 0,
            "logswitch": True,
        }
        line_plot_html = linegraph.plot(self.dumpling_count_plot_data, pconfig)

        # Add a report section with the counts histogram
        self.add_section(
            description="This plot shows the distribution of numbers of observations of variants in a sample.",
            helptext="""
            This plot is useful for determining the quality of a library. In an ideal case,
            all variants will be incorporated into the library at the same frequency. In
            this case, we will expect the observed counts to be poisson distributed with a mean
            equal to the coverage depth. In general, there will be some variability during library
            generation, some selection during library amplification, and some sequencing error.
            These will all contribute to a distribution of counts that will differ from the ideal.

            A major problem for a library would be a large number of missing variants. This can be
            identified here as a spike at 0 counts. We can model this distribution as a zero-inflated
            poisson distribution. We can then use the observed distribution to estimate the parameters
            of the model and determine the probability of missing a variant. This can be used to
            determine whether the library is of sufficient quality to proceed with analysis.
            """,
            plot=line_plot_html,
        )

        # Plot the coverage histogram

        # Create line plot for coverage histogram
        if self.dumpling_coverage_plot_data:
            pconfig = {
                "id": "dumpling_coverage_hist",
                "title": "Coverage histogram",
                "ylab": "Number of positions",
                "xlab": "Coverage (read depth)",
                "ymin": 0,
                "smooth_points": 30,
                "logswitch": True,
            }
            line_plot_html = linegraph.plot(self.dumpling_coverage_plot_data, pconfig)

            # Add a report section for coverage histogram
            self.add_section(
                description="This plot shows the sequencing depth of the target gene.",
                helptext="""
                Help text for coverage histogram.
                """,
                plot=line_plot_html,
            )

        # Reproducibility section is rendered last — it's reference material
        # (versions + config snapshot) rather than QC signal, so it sits at the
        # bottom rather than competing with the actual plots for top placement.
        self._add_reproducibility_section()

    def _add_settings_section(self):
        """Render a table of ORF coordinates and analysis settings from config."""
        rows = []

        if hasattr(config, "orf_start"):
            rows.append(
                f"<tr><td>ORF coordinates</td>"
                f"<td>{config.orf_start}–{config.orf_end} "
                f"({config.orf_length}&nbsp;bp, {config.orf_n_codons}&nbsp;codons)</td></tr>"
            )

        if hasattr(config, "n_variants"):
            rows.append(f"<tr><td>Designed variants</td><td>{config.n_variants}</td></tr>")

        if hasattr(config, "amplicon_length"):
            orf_frac = config.orf_length / config.amplicon_length
            rows.append(f"<tr><td>Amplicon length</td><td>{config.amplicon_length}&nbsp;bp</td></tr>")
            rows.append(
                f"<tr><td>ORF fraction of amplicon</td>"
                f"<td>{orf_frac:.1%} — expected fraction of reads carrying a designed variant</td></tr>"
            )
            rows.append(
                f"<tr><td>Expected WT-like read fraction</td>"
                f"<td>{1.0 - orf_frac:.1%} — amplicon bases outside ORF, always wild-type</td></tr>"
            )

        label_map = {
            "minq": "Minimum base quality (minq)",
            "min_variants": "Minimum variant count filter",
            "amplicon_start": "Amplicon start",
            "amplicon_end": "Amplicon end",
        }
        for key, val in getattr(config, "dumpling_settings", {}).items():
            if key == "amplicon_length":
                continue  # already rendered above
            rows.append(f"<tr><td>{label_map.get(key, key)}</td><td>{val}</td></tr>")

        if not rows:
            return

        html = (
            '<table class="table table-condensed" style="width:auto">'
            "<thead><tr><th>Setting</th><th>Value</th></tr></thead>"
            "<tbody>" + "".join(rows) + "</tbody></table>"
        )

        self.add_section(
            name="Analysis settings",
            anchor="dumpling-settings",
            description="ORF coordinates and analysis parameters for this run.",
            content=html,
        )

    def _add_reproducibility_section(self):
        """Render dumpling version, pipeline-config snapshot, and tool versions.

        Each block is independently optional — we render whichever subset is
        actually available. If none of these are populated, the section is
        skipped entirely.
        """
        rows = []

        version = getattr(config, "dumpling_pipeline_version", None)
        if version:
            rows.append(f"<tr><td>Dumpling pipeline version</td><td><code>{version}</code></td></tr>")

        plugin_version = getattr(config, "multiqc_dumpling_version", None)
        if plugin_version:
            rows.append(f"<tr><td>MultiQC plugin version</td><td><code>{plugin_version}</code></td></tr>")

        # Tool versions discovered while parsing logs.
        for tool in sorted(self.tool_versions):
            rows.append(
                f"<tr><td>{tool}</td><td><code>{self.tool_versions[tool]}</code></td></tr>"
            )

        # Pipeline config snapshot — render a curated subset of high-signal keys.
        pipeline_config = getattr(config, "dumpling_pipeline_config", None)
        if isinstance(pipeline_config, dict):
            reproducibility_keys = [
                "experiment", "orf", "scoring_backend", "enrich2",
                "aligner", "noprocess", "remove_zeros", "regenerate_variants",
                "max_deletion_length", "baseline_condition",
                "min_q", "min_variant_obs", "kmers", "sam",
                "bbtools_compression",
            ]
            for key in reproducibility_keys:
                if key in pipeline_config:
                    val = pipeline_config[key]
                    rows.append(f"<tr><td>config.{key}</td><td><code>{val}</code></td></tr>")

        if not rows:
            return

        html = (
            '<table class="table table-condensed" style="width:auto">'
            "<thead><tr><th>Field</th><th>Value</th></tr></thead>"
            "<tbody>" + "".join(rows) + "</tbody></table>"
        )
        self.add_section(
            name="Reproducibility",
            anchor="dumpling-reproducibility",
            description=(
                "Pipeline and tool versions, plus a snapshot of the dumpling "
                "snakemake config used for this run. Provide "
                "<code>multiqc_dumpling.pipeline_version</code> and "
                "<code>multiqc_dumpling.pipeline_config_file</code> in the "
                "MultiQC config to populate this section."
            ),
            content=html,
        )

    @staticmethod
    def _classify_hgvs(hgvs: str) -> str:
        """Classify a variant by its HGVS protein notation (one-letter amino acid codes).

        Returns one of: wildtype, synonymous, missense, nonsense,
        insertion, deletion, indel, other.

        Dumpling's ``process_variants.name_to_hgvs`` wraps every name in
        parens (``"p.(" + name + ")"``), so real strings look like ``p.(M1F)``.
        The regex tolerates both parenthesized and unparenthesized forms.
        """
        if not isinstance(hgvs, str):
            return "other"
        h = hgvs.strip()
        if h in ("p.=", "p.(=)", "=", "WT", "wt"):
            return "wildtype"
        # Check indel before separate del/ins so "delins" is caught first
        if "delins" in h:
            return "indel"
        if "del" in h:
            return "deletion"
        if "ins" in h:
            return "insertion"
        m = re.match(r"p\.\(?([A-Z*])(\d+)([A-Z*])\)?", h)
        if m:
            orig_aa, new_aa = m.group(1), m.group(3)
            if new_aa == "*":
                return "nonsense"
            if orig_aa == new_aa:
                return "synonymous"
            return "missense"
        return "other"

    def parse_counts(self, file, bin_n=30):
        """Parse the counts file and return a dict of counts and their frequencies for
        plotting a histogram.

        At the same time, set the variants number if it doesn't already exist.
        """
        df = pd.read_csv(file, index_col=False, header=0)

        max_counts = int(df["count"].max())
        mean_counts = float(df["count"].mean())
        median_counts = int(df["count"].median())
        n_zero_counts = int(len(df[df["count"] == 0]))

        if not hasattr(config, "n_variants"):
            log.warning(
                "Number of variants not set: using counts file %s, setting n_variants to %s",
                file,
                len(df),
            )
            config.n_variants = len(df)

        fraction_zero_counts = n_zero_counts / config.n_variants

        counts_stats_dict = {
            "Max counts": max_counts,
            "Mean counts": mean_counts,
            "Median counts": median_counts,
            "Number of zero counts": n_zero_counts,
            "Fraction of zero counts": fraction_zero_counts,
        }

        # Classify designed variants by type using the hgvs column
        if "hgvs" in df.columns:
            type_series = df["hgvs"].apply(self._classify_hgvs)
            type_counts = type_series.value_counts().to_dict()
            total = len(df)
            for vt in _VARIANT_TYPES:
                n = type_counts.get(vt, 0)
                counts_stats_dict[f"n_{vt}"] = n
                counts_stats_dict[f"fraction_{vt}"] = n / total if total > 0 else 0.0

        # If there are no counts, return a dict with a single entry for plotting.

        if max_counts == 0:
            return {0: len(df)}, counts_stats_dict

        # Use 30 bins for the histogram, unless there are fewer than 30 counts, in which case
        # use the number of counts as the number of bins.

        if max_counts < bin_n:
            bin_n = max_counts

        out, bins = pd.cut(
            df["count"], bins=bin_n, include_lowest=True, right=False, retbins=True
        )

        # now make dict between bin and out

        counts_bin_dict = dict(
            zip(bins.tolist(), out.value_counts(normalize=False).sort_index().tolist())
        )

        return counts_bin_dict, counts_stats_dict

    # ──────────────────────────────────────────────────────────────
    # Scoring output parser (rosace / lilace)
    # ──────────────────────────────────────────────────────────────

    @staticmethod
    def _scoring_backend_from_path(root: str) -> str | None:
        """Infer which scoring backend produced a ``*_scores.csv`` from its dir.

        Dumpling writes to ``results/{exp}/rosace/`` or ``results/{exp}/lilace/``;
        we look for those directory names in the path. Returns None when the
        file doesn't look like a scoring output (e.g. a same-named file
        produced by some other tool) — caller should skip it.
        """
        parts = Path(root).parts
        if "rosace" in parts:
            return "rosace"
        if "lilace" in parts:
            return "lilace"
        return None

    def parse_scores(self, file_path) -> tuple[dict, dict, "pd.Series"]:
        """Parse a rosace/lilace per-condition scores CSV.

        Returns ``(stats, histogram, mean_scores)``:
          - **stats**: dict for general-stats columns
              * ``n_total``                — rows in the file
              * ``n_scored``               — variants with non-NA mean AND finite SD
              * ``pct_scored``             — n_scored / n_designed (uses config.n_variants)
              * ``n_significant``          — variants with ``lfsr < LFSR_THRESHOLD``
              * ``fraction_significant``   — n_significant / n_scored
              * ``mean_score``             — mean of the per-variant means
          - **histogram**: dict ``{bin_left: count}`` for the score-distribution plot
          - **mean_scores**: indexed by variant name — used for cross-condition concordance
        """
        df = pd.read_csv(file_path)

        # Required columns. Be loud if missing — that means the format has drifted.
        required = {"variants", "mean", "sd", "lfsr"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Scoring CSV missing required columns: {sorted(missing)}")

        # "Scored" = non-NA mean AND finite SD (user spec). Replace ±inf with NaN
        # before the isna check so values like Inf don't slip through.
        sd_clean = df["sd"].replace([np.inf, -np.inf], np.nan)
        scored_mask = df["mean"].notna() & sd_clean.notna()
        n_scored = int(scored_mask.sum())
        n_total = len(df)

        # Coverage of designed variants — denominator is the designed-variants
        # count from the init hook. If unset, fall back to n_total in the file
        # (the count file's variant set should match) but warn.
        n_designed = getattr(config, "n_variants", None)
        if n_designed is None:
            log.warning(
                "config.n_variants not set — using row count of %s for pct_scored. "
                "Consider configuring multiqc_dumpling.variants_file.",
                file_path,
            )
            n_designed = n_total
        pct_scored = n_scored / n_designed if n_designed else 0.0

        scored = df[scored_mask]
        n_significant = int((scored["lfsr"] < LFSR_THRESHOLD).sum())
        fraction_significant = n_significant / n_scored if n_scored else 0.0

        stats = {
            "n_total": n_total,
            "n_scored": n_scored,
            "pct_scored": pct_scored,
            "n_significant": n_significant,
            "fraction_significant": fraction_significant,
            "mean_score": float(scored["mean"].mean()) if n_scored else 0.0,
        }

        # Histogram of scored means. Use fixed bin count — the score range
        # varies by experiment so we let pandas pick the edges.
        if n_scored:
            counts, bins = np.histogram(scored["mean"].to_numpy(), bins=40)
            histogram = dict(zip(bins[:-1].tolist(), counts.tolist()))
        else:
            histogram = {}

        # Series indexed by variant name — empty index when nothing scored.
        # Used downstream for pairwise concordance between conditions.
        scores = pd.Series(
            scored["mean"].to_numpy(),
            index=scored["variants"].to_numpy(),
            name=str(file_path),
        )

        return stats, histogram, scores

    def parse_enrich2_scores(self, file_path) -> tuple[dict, dict, "pd.Series"]:
        """Parse an enrich2 per-replicate ``main_identifiers_scores.tsv``.

        Enrich2 emits one file per condition/replicate selection
        (``cond_X_R{N}_sel/main_identifiers_scores.tsv``). The file is a flat
        TSV with one row per variant and columns including ``score``, ``SE``,
        and ``pvalue_raw``. We use those three. Returns the same
        ``(stats, histogram, scores)`` triple as :meth:`parse_scores` so the
        downstream report sections don't need to branch on backend.

        Significance is ``pvalue_raw < ENRICH2_PVALUE_THRESHOLD`` (Wald-style
        test of score≠0). This is the enrich2 equivalent of rosace/lilace's
        ``lfsr < LFSR_THRESHOLD`` — same conventional 0.05 cutoff, different
        statistic.
        """
        df = pd.read_csv(file_path, sep="\t")

        required = {"hgvs", "score", "SE", "pvalue_raw"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Enrich2 scores TSV missing required columns: {sorted(missing)}")

        # Same "scored" definition as parse_scores: non-NA score AND finite SE.
        se_clean = df["SE"].replace([np.inf, -np.inf], np.nan)
        scored_mask = df["score"].notna() & se_clean.notna()
        n_scored = int(scored_mask.sum())
        n_total = len(df)

        n_designed = getattr(config, "n_variants", None)
        if n_designed is None:
            log.warning(
                "config.n_variants not set — using row count of %s for pct_scored. "
                "Consider configuring multiqc_dumpling.variants_file.",
                file_path,
            )
            n_designed = n_total
        pct_scored = n_scored / n_designed if n_designed else 0.0

        scored = df[scored_mask]
        n_significant = int((scored["pvalue_raw"] < ENRICH2_PVALUE_THRESHOLD).sum())
        fraction_significant = n_significant / n_scored if n_scored else 0.0

        stats = {
            "n_total": n_total,
            "n_scored": n_scored,
            "pct_scored": pct_scored,
            "n_significant": n_significant,
            "fraction_significant": fraction_significant,
            "mean_score": float(scored["score"].mean()) if n_scored else 0.0,
        }

        if n_scored:
            counts, bins = np.histogram(scored["score"].to_numpy(), bins=40)
            histogram = dict(zip(bins[:-1].tolist(), counts.tolist()))
        else:
            histogram = {}

        scores = pd.Series(
            scored["score"].to_numpy(),
            index=scored["hgvs"].to_numpy(),
            name=str(file_path),
        )

        return stats, histogram, scores

    @staticmethod
    def _enrich2_replicate_name_from_root(root: str) -> str | None:
        """Derive a replicate sample name from the enrich2 ``_sel`` directory.

        ``results/{exp}/enrich/tsv/cond_A_R1_sel/main_identifiers_scores.tsv``
        gives ``cond_A_R1`` (parent dir minus the ``_sel`` suffix). Returns
        None if the parent dir doesn't end in ``_sel`` — guards against the
        ``_exp`` directory or future enrich2 layouts.
        """
        parent = Path(root).name
        if parent.endswith("_sel"):
            return parent[: -len("_sel")]
        return None

    # ──────────────────────────────────────────────────────────────
    # Pipeline step log parsers
    # ──────────────────────────────────────────────────────────────

    @staticmethod
    def _parse_bbduk_log(content: str) -> dict:
        """Parse a BBDuk stderr log (identical format for decontam and trim steps).

        Returns keys: reads_in, bases_in, reads_out, bases_out.
        """
        stats = {}
        m = re.search(r"Input:\s+([\d,]+)\s+reads\s+([\d,]+)\s+bases", content)
        if m:
            stats["reads_in"] = int(m.group(1).replace(",", ""))
            stats["bases_in"] = int(m.group(2).replace(",", ""))
        m = re.search(r"Result:\s+([\d,]+)\s+reads.*?([\d,]+)\s+bases", content)
        if m:
            stats["reads_out"] = int(m.group(1).replace(",", ""))
            stats["bases_out"] = int(m.group(2).replace(",", ""))
        return stats

    @staticmethod
    def _parse_bbmerge_log(content: str) -> dict:
        """Parse a BBMerge stderr log.

        Returns keys:
          pairs_in  — number of read pairs submitted to BBMerge
          reads_in  — pairs_in × 2 (raw read count, for consistency with other parsers)
          reads_out — number of successfully merged reads
          pct_joined

        BBMerge has two output styles depending on version:
          "Pairs:  N"  then  "Joined:  M  P%"
          "Input:  2N reads"  then  "Joined:  M reads (P%)"
        The Joined regex handles both by skipping non-numeric text between count and %.
        """
        stats = {}
        m = re.search(r"Pairs:\s+([\d,]+)", content)
        if m:
            stats["pairs_in"] = int(m.group(1).replace(",", ""))
            stats["reads_in"] = stats["pairs_in"] * 2
        else:
            m = re.search(r"Input:\s+([\d,]+)\s+reads", content)
            if m:
                stats["reads_in"] = int(m.group(1).replace(",", ""))
                stats["pairs_in"] = stats["reads_in"] // 2
        # Match "Joined:  N  P%" or "Joined:  N reads (P%)" — .*? skips any intervening text
        m = re.search(r"Joined:\s+([\d,]+).*?(\d+\.?\d+)%", content)
        if m:
            stats["reads_out"] = int(m.group(1).replace(",", ""))
            stats["pct_joined"] = float(m.group(2))
        return stats

    # ──────────────────────────────────────────────────────────────
    # Tool version extraction (reproducibility)
    # ──────────────────────────────────────────────────────────────

    # Regex map: tool label → regex to find its version string in a log file.
    # Patterns return the version in group(1).
    _TOOL_VERSION_PATTERNS = {
        "BBTools":  re.compile(r"^Version\s+(\S+)", re.MULTILINE),
        "GATK":     re.compile(r"Genome Analysis Toolkit \(GATK\) v(\S+)"),
        "samtools": re.compile(r"produced by samtools \S+ \(([\d.+a-z\-]+)\)"),
    }

    @classmethod
    def _extract_tool_version(cls, content: str, tool: str) -> str | None:
        """Pull the first version string for ``tool`` out of a log's content.

        Returns None if not found. Used to populate the reproducibility section
        without re-reading files — we run it on logs already loaded for parsing.
        """
        pat = cls._TOOL_VERSION_PATTERNS.get(tool)
        if pat is None:
            return None
        m = pat.search(content)
        return m.group(1) if m else None

    @staticmethod
    def _parse_samtools_stats(content: str) -> dict:
        """Parse a ``samtools stats`` output for the mapping-step summary only.

        Extracts reads_in/reads_out for the read waterfall and bases_in/bases_out
        for the bases waterfall. ``total length`` is the input base count;
        ``bases mapped (cigar)`` is the more accurate of two "bases mapped"
        values samtools reports (the comment in the file: ``# more accurate``).
        MultiQC's built-in samtools module already renders the full section in
        the report; this parser is just for the waterfall plots.
        """
        stats = {}
        wanted = {
            "raw total sequences": "reads_in",
            "reads mapped":        "reads_out",
            "total length":        "bases_in",
            "bases mapped (cigar)": "bases_out",
        }
        for line in content.splitlines():
            if not line.startswith("SN\t"):
                continue
            # Lines look like:  SN<TAB>key:<TAB>value<TAB>#comment
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].rstrip(":").strip()
            if key not in wanted:
                continue
            try:
                stats[wanted[key]] = int(parts[2].strip())
            except ValueError:
                continue
        if "reads_in" in stats and stats["reads_in"] and "reads_out" in stats:
            stats["pct_mapped"] = stats["reads_out"] / stats["reads_in"] * 100
        return stats

    @staticmethod
    def _parse_bbmap_log(content: str) -> dict:
        """Parse a BBMap stderr log.

        Returns keys: reads_in, reads_out (mapped), pct_mapped, plus
        bases_in, bases_out when those lines are present.

        BBMap has two summary styles depending on flags/version:
          1. ``Mapped:  N  P%``  (single-line, used in early test fixtures)
          2. ``Mapped reads:  N`` / ``Mapped bases:  N`` (the form dumpling
             actually emits with covstats/etc. requested).
        Both are matched here.
        """
        stats = {}
        # "Reads Used:  N	(M bases)" gives reads_in and bases_in in one line
        m = re.search(r"Reads Used:\s+([\d,]+)\s*\(?([\d,]+)?\s*bases\)?", content)
        if m:
            stats["reads_in"] = int(m.group(1).replace(",", ""))
            if m.group(2):
                stats["bases_in"] = int(m.group(2).replace(",", ""))

        m = re.search(r"^Mapped:\s+([\d,]+)\s+([\d.]+)%", content, re.MULTILINE)
        if m:
            stats["reads_out"] = int(m.group(1).replace(",", ""))
            stats["pct_mapped"] = float(m.group(2))
            return stats

        m = re.search(r"^Mapped reads:\s+([\d,]+)", content, re.MULTILINE)
        if m:
            stats["reads_out"] = int(m.group(1).replace(",", ""))
        m = re.search(r"^Mapped bases:\s+([\d,]+)", content, re.MULTILINE)
        if m:
            stats["bases_out"] = int(m.group(1).replace(",", ""))
        m = re.search(r"^Percent mapped:\s+([\d.]+)", content, re.MULTILINE)
        if m:
            stats["pct_mapped"] = float(m.group(1))
        return stats

    @staticmethod
    def _parse_kv_stats_file(content: str) -> dict:
        """Parse a dumpling processing-stats TSV (total / accepted / rejected).

        These files are written by ``process_variants.write_stats_file`` as
        long-format ``key\\tvalue`` lines with no header row. Values are
        integers; non-integer values are skipped with a debug log.

        The same shape is used for ``*_total_processing.tsv``,
        ``*_accepted_processing.tsv`` and ``*_rejected_processing.tsv``, so
        callers don't need different parsers per file type.
        """
        stats = {}
        for raw_line in content.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                log.debug("Skipping malformed stats line: %r", raw_line)
                continue
            key, value = parts[0].strip(), parts[1].strip()
            try:
                stats[key] = int(value)
            except ValueError:
                log.debug("Non-integer stats value for %s: %r", key, value)
        return stats

    def _compute_processing_metrics(self) -> dict:
        """Derive per-sample ratios from the total/accepted/rejected stats dicts.

        Returns one dict per sample keyed by metric name. Skips samples that
        don't have a ``total_counts`` value (the denominator), and logs a
        warning when ``total_accepted + total_rejected != total_counts`` —
        which would indicate a pipeline accounting bug.
        """
        all_samples = set(self.pipeline_filter_data) | set(self.pipeline_rejected_data)
        out = {}
        for sample in all_samples:
            total = self.pipeline_filter_data.get(sample, {})
            rejected = self.pipeline_rejected_data.get(sample, {})

            total_counts = total.get("total_counts")
            if not total_counts:
                continue

            accepted = total.get("total_accepted_counts", 0)
            rejected_total = total.get("total_rejected_counts", 0)

            # Sanity check: pipeline invariant accepted + rejected == total.
            # If broken, surface it loudly rather than letting downstream ratios silently
            # add to >100%.
            if accepted + rejected_total != total_counts:
                log.warning(
                    "Sample %s: accepted (%d) + rejected (%d) != total (%d) — "
                    "pipeline stats inconsistent",
                    sample, accepted, rejected_total, total_counts,
                )

            metrics = {
                "read_acceptance_rate": accepted / total_counts,
                "pct_reads_outside_orf": rejected.get("outside_orf_counts", 0) / total_counts,
                "pct_reads_frameshift": rejected.get("fs_counts", 0) / total_counts,
                "pct_reads_wrong_variant": rejected.get("wrong_variant_counts", 0) / total_counts,
            }
            out[sample] = metrics
        return out

    # ──────────────────────────────────────────────────────────────
    # Scoring report sections
    # ──────────────────────────────────────────────────────────────

    def _add_scoring_general_stats(self):
        """Promote scoring summary metrics into the general-stats table."""
        headers = {
            "pct_scored": {
                "title": "% scored",
                "description": (
                    f"Fraction of designed variants with a non-NA mean and finite SD "
                    f"(threshold lfsr<{LFSR_THRESHOLD} not applied here)."
                ),
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
                "format": "{:.1%}",
            },
            "fraction_significant": {
                "title": "% sig",
                "description": (
                    f"Fraction of scored variants below the conventional 0.05 significance "
                    f"cutoff. Rosace/lilace use lfsr<{LFSR_THRESHOLD}; enrich2 uses "
                    f"pvalue_raw<{ENRICH2_PVALUE_THRESHOLD} (Wald-style test of score≠0)."
                ),
                "min": 0,
                "max": 1,
                "scale": "Blues",
                "format": "{:.1%}",
            },
            "n_significant": {
                "title": "# significant",
                "description": "Number of significant variants at the conventional 0.05 cutoff (per-backend metric — see % sig).",
                "min": 0,
                "scale": "Blues",
                "format": "{:,}",
                "hidden": True,
            },
            "mean_score": {
                "title": "Mean score",
                "description": "Mean of per-variant score (rosace/lilace ``mean`` column).",
                "scale": "RdBu",
                "format": "{:.2f}",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.scoring_stats, headers)

    def _add_scoring_histogram_section(self):
        """Per-condition score-distribution histogram (one line per condition)."""
        plot_data = {k: v for k, v in self.scoring_plot_data.items() if v}
        if not plot_data:
            return
        pconfig = {
            "id": "dumpling_scoring_histogram",
            "title": "Variant score distributions",
            "ylab": "Number of variants",
            "xlab": "Mean score (rosace/lilace)",
            "ymin": 0,
        }
        self.add_section(
            name="Variant scores",
            anchor="dumpling-scoring-histogram",
            description=(
                "Distribution of per-variant mean scores from rosace/lilace. "
                "In a typical selection screen, synonymous variants cluster near zero "
                "(neutral) and missense variants spread to either side."
            ),
            plot=linegraph.plot(plot_data, pconfig),
        )

    def _add_scoring_concordance_section(self):
        """Cross-condition concordance — pairwise Pearson + Spearman between conditions.

        With multiple conditions, replicates of similar phenotypes should
        correlate more strongly than dissimilar ones; this makes
        within-vs-between comparisons explicit.
        """
        keys = list(self.scoring_score_vectors.keys())
        if len(keys) < 2:
            return  # Need at least 2 score files to compare

        # Align scores by variant name. ``pd.concat`` does an outer join on
        # the index; per-pair correlation drops NaNs (the union of variants
        # is the right baseline for "how much do these two agree").
        aligned = pd.concat(self.scoring_score_vectors, axis=1, join="outer")

        rows = []
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                pair = aligned[[a, b]].dropna()
                if len(pair) < 3:
                    continue  # Not enough overlap to compute a meaningful correlation
                pearson = pair[a].corr(pair[b], method="pearson")
                spearman = pair[a].corr(pair[b], method="spearman")
                rows.append((a, b, len(pair), pearson, spearman))

        if not rows:
            return

        body = "".join(
            f"<tr><td>{a}</td><td>{b}</td><td>{n:,}</td>"
            f"<td>{p:.3f}</td><td>{s:.3f}</td></tr>"
            for a, b, n, p, s in rows
        )
        html = (
            '<table class="table table-condensed" style="width:auto">'
            "<thead><tr><th>Condition A</th><th>Condition B</th>"
            "<th>n shared</th><th>Pearson r</th><th>Spearman ρ</th></tr></thead>"
            f"<tbody>{body}</tbody></table>"
        )
        self.add_section(
            name="Score concordance",
            anchor="dumpling-scoring-concordance",
            description=(
                "Pairwise correlation of per-variant mean scores between scoring outputs. "
                "Computed over variants scored in both files (NaNs dropped pairwise). "
                "High concordance between phenotypically similar conditions is expected; "
                "low concordance flags possible technical issues."
            ),
            content=html,
        )

    def _add_pipeline_waterfall_section(self):
        """Render the combined read+base pipeline survival waterfall.

        One linegraph with two switchable datasets (reads in pair-equivalent
        fragments; total bases). Step order mirrors the actual pipeline:
        trim → decontam → EC → map → filter. Per-step unit conversions are
        documented inline in ``read_specs`` and ``base_specs``. The summary
        tables sit behind ``<details>`` toggles, closed by default, so the
        section reads as just plot + brief description.

        TODO: validate on data with realistic adapter content. The current
        example_experiment fixture trims <0.5% of bases, so the bases plot
        looks almost flat. A library with proper adapter contamination will
        show a much more dramatic step at the trim stage.
        """
        # (label, source dict, key, conversion to fragment-equivalents)
        #
        # Unit notes — every value here is in *fragments* (a pair counts as 1):
        #   * BBDuk Input/Result lines report ``reads`` (mates), so ÷2 for pairs.
        #   * BBMerge in dumpling runs as ``ecco mix`` — error-correct only,
        #     pairs stay as pairs. Fragment count is unchanged, so we use
        #     ``pairs_in`` (not the ``Joined`` count, which is the subset that
        #     *could* have been merged but in ecco mode are kept as pairs).
        #   * BBMap "Reads Used"/"Mapped reads" and samtools "raw total
        #     sequences"/"reads mapped" count paired mates separately
        #     (1 pair = 2). ÷2 to get fragment-equivalents.
        #   * The variant-filter step counts variant observations, which in DMS
        #     is ~1 per mapped fragment.
        read_specs = [
            ("Input",     self.pipeline_trim_data,     "reads_in",              lambda v: v / 2),
            ("Trimmed",   self.pipeline_trim_data,     "reads_out",             lambda v: v / 2),
            ("Decontam",  self.pipeline_decontam_data, "reads_out",             lambda v: v / 2),
            ("EC'd",      self.pipeline_merge_data,    "pairs_in",              lambda v: v),
            ("Mapped",    self.pipeline_map_data,      "reads_out",             lambda v: v / 2),
            ("Var-pass",  self.pipeline_filter_data,   "total_accepted_counts", lambda v: v),
        ]
        # Bases: no ÷2 (both mates contribute). EC step reuses decontam output
        # because the bbmerge log has no base count and ``ecco`` doesn't drop
        # bases. Variant-filter is omitted — dumpling tracks variant counts at
        # that stage, not bases.
        base_specs = [
            ("Input",    self.pipeline_trim_data,     "bases_in",  lambda v: v),
            ("Trimmed",  self.pipeline_trim_data,     "bases_out", lambda v: v),
            ("Decontam", self.pipeline_decontam_data, "bases_out", lambda v: v),
            ("EC'd",     self.pipeline_decontam_data, "bases_out", lambda v: v),
            ("Mapped",   self.pipeline_map_data,      "bases_out", lambda v: v),
        ]

        read_data = self._build_waterfall_series(read_specs)
        base_data = self._build_waterfall_series(base_specs)
        if not read_data and not base_data:
            return

        pconfig = {
            "id": "dumpling_pipeline_waterfall",
            "title": "Dumpling: Pipeline survival",
            "xlab": "Pipeline step",
            "ymin": 0,
            "logswitch": True,
            "categories": True,
            "data_labels": [
                {"name": "Reads (fragments)", "ylab": "Fragments (pair-equivalent)"},
                {"name": "Bases",             "ylab": "Bases"},
            ],
        }
        plot_html = linegraph.plot([read_data, base_data], pconfig)

        # Plot ↔ Table toggle — MultiQC linegraph has no native "Table" view
        # (unlike bargraph's cpswitch), so we render our own button group +
        # hidden table div and toggle visibility with a small inline script.
        # Styling uses Bootstrap classes already loaded by MultiQC so it
        # matches the look of the existing Log10 / dataset-switcher buttons.
        read_table = self._build_waterfall_table(
            read_data, [label for label, *_ in read_specs]
        )
        base_table = self._build_waterfall_table(
            base_data, [label for label, *_ in base_specs]
        )
        content_before_plot = (
            '<div class="dumpling-waterfall-toggle" style="text-align:right; '
            'margin-bottom:0.5em">'
            '<div class="btn-group btn-group-sm" role="group">'
            '<button type="button" class="btn btn-outline-secondary active" '
            'data-mode="plot">Plot</button>'
            '<button type="button" class="btn btn-outline-secondary" '
            'data-mode="table">Table</button>'
            "</div></div>"
        )
        content = (
            '<div id="dumpling-waterfall-tables" style="display:none; margin-top:1em">'
            "<h5 style=\"margin-top:0.5em\">Reads (fragments)</h5>"
            f"{read_table}"
            '<h5 style="margin-top:1em">Bases</h5>'
            f"{base_table}"
            "</div>"
            # Inline toggle script — vanilla JS, scoped to this section.
            "<script>(function(){"
            "var sec=document.getElementById('mqc-section-wrapper-dumpling-pipeline-waterfall');"
            "if(!sec)return;"
            "var plot=sec.querySelector('#dumpling_pipeline_waterfall-wrapper');"
            "var tables=sec.querySelector('#dumpling-waterfall-tables');"
            "var btns=sec.querySelectorAll('.dumpling-waterfall-toggle button');"
            "btns.forEach(function(btn){btn.addEventListener('click',function(){"
            "btns.forEach(function(b){b.classList.remove('active');});"
            "btn.classList.add('active');"
            "var mode=btn.getAttribute('data-mode');"
            "if(mode==='plot'){if(plot)plot.style.display='';tables.style.display='none';}"
            "else{if(plot)plot.style.display='none';tables.style.display='';}"
            "});});"
            "})();</script>"
        )

        # autoformat=False to preserve the inline <script> tag — MultiQC's
        # markdown autoformat sanitizes scripts out otherwise. The description
        # is pre-rendered HTML to compensate.
        description_html = (
            "Reads (fragments) and bases remaining at each pipeline step. "
            "Use the <strong>Plot / Table</strong> toggle (top-right) to switch "
            "between the plot and the per-step survival tables; within the plot, "
            "the <strong>Reads (fragments) / Bases</strong> switch picks the dataset. "
            "Reads are in pair-equivalent units (1 pair = 1); bases are summed across "
            "both mates. The EC step (BBMerge <code>ecco mix</code>) only "
            "error-corrects pairs without joining them, so fragment and base counts "
            "are unchanged there by design."
        )
        self.add_section(
            name="Pipeline survival",
            anchor="dumpling-pipeline-waterfall",
            description=description_html,
            content_before_plot=content_before_plot,
            content=content,
            plot=plot_html,
            autoformat=False,
        )

    def _build_waterfall_series(self, step_specs):
        """Build ``{sample: {step_label: value}}`` from a step_specs list.

        Returns an empty dict if no source has any samples.
        """
        all_samples = set()
        for _label, source, _key, _convert in step_specs:
            all_samples.update(source.keys())
        if not all_samples:
            return {}

        out: dict[str, dict[str, float]] = {}
        for sample in sorted(all_samples):
            series: dict[str, float] = {}
            for label, source, key, convert in step_specs:
                raw = source.get(sample, {}).get(key)
                if raw is None:
                    continue
                series[label] = convert(raw)
            if series:
                out[sample] = series
        return out

    @staticmethod
    def _build_waterfall_table(series_by_sample, step_labels):
        """Render the absolute-value + step-relative-survival HTML table."""
        rows_html = []
        for sample, series in series_by_sample.items():
            cells = [f"<td><b>{sample}</b></td>"]
            prev = None
            for label in step_labels:
                val = series.get(label)
                if val is None:
                    cells.append("<td>—</td>")
                    continue
                txt = f"{val:,.0f}"
                if prev:
                    survival = val / prev * 100
                    color = "#c0392b" if survival < 80 else "inherit"
                    txt += f'<br><small style="color:{color}">{survival:.1f}%</small>'
                cells.append(f"<td>{txt}</td>")
                prev = val
            rows_html.append("<tr>" + "".join(cells) + "</tr>")
        th = "".join(f"<th>{s}</th>" for s in step_labels)
        return (
            '<div style="overflow-x:auto"><table class="table table-condensed table-bordered">'
            f"<thead><tr><th>Sample</th>{th}</tr></thead>"
            f"<tbody>{''.join(rows_html)}</tbody>"
            "</table></div>"
        )

    def parse_coverage(self, file, bin_n=300):
        """Parse the coverage file and return a dict of coverage and their frequencies for
        plotting a histogram. Coverage is calculated by scaling the observed coverage by the
        number of variants.

        The coverage file is a tab delimited file with two columns: position in orf (RefPos)
        and coverage.
        """
        df = pd.read_csv(file, index_col=False, header=0, sep="\t")

        # Calculate the coverage by scaling the observed coverage by the number of variants.
        # The number of variants is the number of rows in the counts file.

        df["Coverage"] = df["Coverage"] / config.n_variants
        df["Coverage"] = df["Coverage"].round(0).astype(int)

        # Use 30 bins for the histogram, unless there are fewer than 30 counts, in which case
        # use the number of counts as the number of bins.

        max_coverage = int(df["Coverage"].max())
        min_coverage = int(df["Coverage"].min())
        mean_coverage = float(df["Coverage"].mean())

        bin_n = min(bin_n, max_coverage)
        bin_n = max(1, int(bin_n))

        logging.debug("Using %s bins", bin_n)

        out, bins = pd.cut(
            df["Coverage"], bins=bin_n, include_lowest=True, right=False, retbins=True
        )

        # now make dict between bin and out

        coverage_bin_dict = dict(
            zip(bins.tolist(), out.value_counts(normalize=False).sort_index().tolist())
        )

        coverage_stats_dict = {
            "Max coverage": max_coverage,
            "Min coverage": min_coverage,
            "Mean coverage": mean_coverage,
        }

        return coverage_bin_dict, coverage_stats_dict
