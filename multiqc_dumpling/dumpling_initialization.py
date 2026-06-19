#!/usr/bin/env python

from __future__ import print_function
import importlib_metadata
from pathlib import Path

import logging

from multiqc import config

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

__version__ = importlib_metadata.version("multiqc_dumpling")
config.multiqc_dumpling_version = __version__


def _populate_fastq_sample_renames(experiment_path: Path) -> None:
    """Read dumpling's experiment CSV and register fastq→sample rename rules.

    The CSV is expected to have ``sample`` and ``file`` columns, where
    ``file`` is the FASTQ basename minus the ``_R[12]_001.fastq.gz`` suffix.
    We register substring replacements so any sample name containing the
    FASTQ prefix (e.g. ``1_S1_L001_R1_001``) gets the prefix replaced with
    the dumpling sample name (``A_R1_T0``); MultiQC's existing fastqc clean
    extensions then strip the residual ``_R1_001`` suffix.

    Ambiguous fastq prefixes — those used by more than one sample row in
    the CSV — are deliberately *not* renamed. Renaming such a prefix would
    arbitrarily credit one sample with QC data that belongs to several, so
    we leave those rows as-is and log a warning instead.
    """
    import csv

    # First pass: build a multimap fastq_prefix → [sample, ...] so we can
    # detect M:N cases before deciding what to rename.
    prefix_to_samples: dict[str, list[str]] = {}
    with experiment_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "sample" not in reader.fieldnames or "file" not in reader.fieldnames:
            log.warning(
                "experiment_file %s missing required 'sample' or 'file' column; "
                "skipping fastq→sample rename setup.",
                experiment_path,
            )
            return
        for row in reader:
            sample = (row.get("sample") or "").strip()
            fastq_prefix = (row.get("file") or "").strip()
            if not sample or not fastq_prefix:
                continue
            if sample == fastq_prefix:
                continue  # identity rename — nothing to do
            prefix_to_samples.setdefault(fastq_prefix, []).append(sample)

    rename_pairs: dict[str, str] = {}
    ambiguous: dict[str, list[str]] = {}
    for prefix, samples in prefix_to_samples.items():
        if len(set(samples)) == 1:
            rename_pairs[prefix] = samples[0]
        else:
            ambiguous[prefix] = sorted(set(samples))

    if ambiguous:
        log.warning(
            "Skipping fastq→sample rename for %d FASTQ prefix(es) shared by multiple "
            "samples (per-fastq QC rows will be kept as-is to avoid mis-attribution): %s",
            len(ambiguous),
            ", ".join(f"{p} → {{{', '.join(s)}}}" for p, s in ambiguous.items()),
        )

    if not rename_pairs:
        return

    existing = dict(getattr(config, "sample_names_replace", {}) or {})
    for pattern, replacement in rename_pairs.items():
        existing.setdefault(pattern, replacement)
    config.sample_names_replace = existing
    log.info(
        "Registered %d fastq→sample rename rules from %s",
        len(rename_pairs), experiment_path,
    )


def dumpling_plugin_execution_start():
    """Code to initialize the dumpling MultiQC plugin.

    This adds search patterns for counts outputs and parses
    the snakemake config file to set some parameters of the target.
    """
    log.info("Running dumpling MultiQC Plugin v%s", config.multiqc_dumpling_version)

    # Add to the search patterns used by modules

    if "dumpling/counts" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "dumpling/counts": {
                    "fn": "*.csv",
                    "contents_re": ".*mutation,length,hgvs.*",
                    "num_lines": 2,
                }
            },
        )
    if "gatk/analyze_saturation_mutagenesis/refcoverage" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "gatk/analyze_saturation_mutagenesis/refcoverage": {
                    "fn": "*.refCoverage"
                }
            },
        )
    if "gatk/analyze_saturation_mutagenesis/coveragelengthcounts" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "gatk/analyze_saturation_mutagenesis/coveragelengthcounts": {
                    "fn": "*.coverageLengthCounts"
                }
            },
        )
    # File patterns accept both the legacy snake-case names (kept for any
    # older runs still around) and the current dumpling naming. MultiQC's
    # SearchPattern only accepts a single string for ``fn``, so we use
    # ``fn_re`` (regex) when we need to match more than one filename style.
    if "dumpling/bbduk_decontam" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/bbduk_decontam": {"fn_re": r".+(?:\.clean\.bbduk|_decontam)\.log$"}},
        )
    if "dumpling/bbduk_trim" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/bbduk_trim": {"fn_re": r".+(?:\.trim\.bbduk|_trim)\.log$"}},
        )
    if "dumpling/bbmerge" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/bbmerge": {"fn_re": r".+(?:\.bbmerge|_merge)\.log$"}},
        )
    if "dumpling/bbmap" not in config.sp:
        # ``*_map.log`` matches both ``baseline_1_map.log`` (legacy) and
        # ``baseline_1.bbmap_map.log`` (current) because ``*`` is greedy.
        config.update_dict(
            config.sp,
            {"dumpling/bbmap": {"fn": "*_map.log"}},
        )
    # samtools stats: produced by both aligner branches in the current
    # pipeline (bbmap runs it on the resulting BAM; minimap2 uses it as
    # the primary QC output). Preferring this lets the waterfall work
    # uniformly regardless of which aligner the run used.
    if "dumpling/samtools_stats" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/samtools_stats": {"fn": "*_samtools_stats.txt"}},
        )

    # Strip dumpling-specific filename suffixes so each tool's output for a
    # given sample collapses to the same s_name. Without this, samtools stats
    # for ``baseline_1`` would land under ``baseline_1_samtools_stats`` and
    # the waterfall wouldn't merge with the bbduk/bbmerge data.
    # Order matters — MultiQC strips the first matching suffix. We list the
    # longer/more-specific exts first so e.g. ``.bbmap_map`` is stripped before
    # a shorter ``_map`` (which a user might have added to their own config)
    # gets a chance to take a bite and leave ``.bbmap`` as a leftover suffix.
    # ``.bbmap`` is included as a safety net in case the order doesn't hold.
    _dumpling_clean_exts = [
        ".bbmap_map",
        ".clean.bbduk",
        ".trim.bbduk",
        ".bbmerge",
        ".bbmap",
        ".refCoverage",
        "_samtools_stats",
        "_samtools_flagstat",
        "_total_processing",
        "_accepted_processing",
        "_rejected_processing",
        "_scores",
        # FASTQ read-direction suffixes ("_R1_001" / "_R2_001") aren't in
        # MultiQC's defaults but they're how the Illumina pipeline names
        # files. Stripping them BEFORE sample_names_replace runs lets the
        # cleaned name (the bare fastq prefix) match a rename rule directly,
        # otherwise the rename produces a name like ``A_R1_T3_R1_001``.
        "_R1_001",
        "_R2_001",
    ]
    existing = list(getattr(config, "fn_clean_exts", []) or [])
    for ext in _dumpling_clean_exts:
        if ext not in existing:
            existing.append(ext)
    config.fn_clean_exts = existing
    # ``shared: True`` lets these files be picked up by the plugin AND by any
    # user-defined ``custom_data`` patterns that also match the same filenames
    # — without it, whichever search-key was registered first claims the file
    # and the other gets nothing. The dumpling-side multiqc_config.yaml ships
    # such custom_data entries for the bargraph views, so this is required for
    # both views to co-exist.
    if "dumpling/filter_stats" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/filter_stats": {"fn": "*_total_processing.tsv", "shared": True}},
        )
    if "dumpling/accepted_stats" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/accepted_stats": {"fn": "*_accepted_processing.tsv", "shared": True}},
        )
    if "dumpling/rejected_stats" not in config.sp:
        config.update_dict(
            config.sp,
            {"dumpling/rejected_stats": {"fn": "*_rejected_processing.tsv", "shared": True}},
        )

    # The dumpling-side multiqc_config.yaml defines custom_data entries with
    # these search keys that match the same processing-stats filenames.
    # MultiQC normally lets the first matching non-shared pattern claim the
    # file exclusively; adding those keys to ``filesearch_file_shared`` makes
    # the file flow through to the plugin as well, so both the bargraph view
    # and the plugin's derived-metric view co-exist.
    _shared_keys = ("variant_processing_total", "variant_processing_accepted",
                    "variant_processing_rejected")
    existing_shared = list(getattr(config, "filesearch_file_shared", []) or [])
    for key in _shared_keys:
        if key not in existing_shared:
            existing_shared.append(key)
    config.filesearch_file_shared = existing_shared
    # Scoring outputs: rosace and lilace emit the same per-condition CSV
    # (one row per variant with mean/sd/lfsr columns). The header is the
    # uniquely-identifying string we match on; the plugin distinguishes
    # rosace vs lilace at parse time by the parent directory name.
    if "dumpling/scores" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "dumpling/scores": {
                    "fn": "*_scores.csv",
                    "contents_re": r'^"?variants"?,"?position"?,.*"?mean"?,"?sd"?,"?lfsr"?',
                    "num_lines": 1,
                }
            },
        )
    # Enrich2 per-replicate scoring: dumpling emits one
    # ``{condition}_R{N}_sel/main_identifiers_scores.tsv`` per replicate.
    # The combined experiment-level outputs (``_shared.tsv``, ``_pvalues_wt.tsv``)
    # use a different multi-level-header format and aren't parsed here yet —
    # per-replicate files give us the same information plus replicate
    # concordance, so this is the better starting point.
    if "dumpling/enrich2_scores" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "dumpling/enrich2_scores": {
                    "fn": "main_identifiers_scores.tsv",
                    "contents_re": r"^hgvs\tscore\tSE\b",
                    "num_lines": 1,
                }
            },
        )

    plugin_config = getattr(config, "multiqc_dumpling", None)
    if not isinstance(plugin_config, dict):
        # No dumpling config section — this is not a dumpling run; skip library-specific setup.
        log.debug("No 'multiqc_dumpling' config section found; skipping dumpling-specific initialization.")
        return
    if "orf" not in plugin_config:
        raise ValueError(
            "Missing required 'multiqc_dumpling.orf' setting in MultiQC config."
        )
    if "variants_file" not in plugin_config:
        raise ValueError(
            "Missing required 'multiqc_dumpling.variants_file' setting in MultiQC config."
        )

    # Parse ORF coordinates — expected format: "XX-YY" (1-based, inclusive)
    coords = plugin_config["orf"].split("-")
    if len(coords) != 2:
        log.error("Invalid ORF coordinates: %s", plugin_config["orf"])
        raise ValueError("Invalid ORF coordinates")
    if not all(coord.isdigit() for coord in coords):
        log.error("Invalid ORF coordinates: %s", plugin_config["orf"])
        raise ValueError("Invalid ORF coordinates")
    if int(coords[0]) > int(coords[1]):
        log.error("Invalid ORF coordinates: %s", plugin_config["orf"])
        raise ValueError("Invalid ORF coordinates")
    if int(coords[0]) < 1:
        log.error("Invalid ORF coordinates: %s", plugin_config["orf"])
        raise ValueError("Invalid ORF coordinates")
    start = int(coords[0])
    end = int(coords[1])
    config.orf_start = start
    config.orf_end = end
    config.orf_length = end - start + 1
    config.orf_n_codons = config.orf_length // 3

    # Read optional analysis settings from config
    analysis_settings = {}
    int_keys = ["amplicon_length", "amplicon_start", "amplicon_end", "minq", "min_variants"]
    for key in int_keys:
        if key in plugin_config:
            try:
                val = int(plugin_config[key])
                analysis_settings[key] = val
                setattr(config, key, val)
            except (ValueError, TypeError):
                log.warning("Invalid value for 'multiqc_dumpling.%s': %s", key, plugin_config[key])
    config.dumpling_settings = analysis_settings

    # Reproducibility metadata. Both keys are optional; when present, the
    # module renders them in the Analysis settings section.
    #   pipeline_version       — string like "v0.2.0-5-ge6c2862" (git describe output)
    #   pipeline_config_file   — path to the snakemake config used for the run
    config.dumpling_pipeline_version = plugin_config.get("pipeline_version")
    pipeline_config_file = plugin_config.get("pipeline_config_file")
    config.dumpling_pipeline_config = None
    if pipeline_config_file:
        cfg_path = Path(pipeline_config_file)
        if cfg_path.is_file():
            try:
                import yaml
                with cfg_path.open("r") as f:
                    config.dumpling_pipeline_config = yaml.safe_load(f)
            except Exception as exc:
                log.warning("Could not parse multiqc_dumpling.pipeline_config_file %s: %s", cfg_path, exc)
        else:
            log.warning("multiqc_dumpling.pipeline_config_file does not exist: %s", cfg_path)

    # Experiment-file plumbing — fastqc names its output after the FASTQ
    # filename prefix (e.g. ``1_S1_L001_R1_001``), but every downstream step
    # uses dumpling's sample name (e.g. ``A_R1_T0``). The experiment CSV
    # records the mapping in its ``file → sample`` columns. We use those to
    # populate ``config.sample_names_replace`` so MultiQC merges fastqc rows
    # into the same general-stats row as the rest of the per-sample data.
    experiment_file = plugin_config.get("experiment_file")
    if experiment_file:
        exp_path = Path(experiment_file)
        if exp_path.is_file():
            _populate_fastq_sample_renames(exp_path)
        else:
            log.warning("multiqc_dumpling.experiment_file does not exist: %s", exp_path)

    # Read the designed variants file to get the number of designed variants
    # The first line is a header, so subtract 1 from the total number of lines

    variants_path = Path(plugin_config["variants_file"])
    if not variants_path.is_file():
        raise FileNotFoundError(
            "Configured multiqc_dumpling.variants_file does not exist: "
            f"{plugin_config['variants_file']}"
        )

    with variants_path.open("r") as f:
        config.n_variants = sum(1 for line in f) - 1

    # Add the designed variants file to the ignore list
    config.fn_ignore_files.extend([str(variants_path)])
