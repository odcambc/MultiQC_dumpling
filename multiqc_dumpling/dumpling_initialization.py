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


def dumpling_plugin_execution_start():
    """Code to intialize the dumpling MultiQC plugin.

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

    plugin_config = getattr(config, "multiqc_dumpling", None)
    if not isinstance(plugin_config, dict):
        raise ValueError(
            "Missing required 'multiqc_dumpling' configuration section in MultiQC config."
        )
    if "orf" not in plugin_config:
        raise ValueError(
            "Missing required 'multiqc_dumpling.orf' setting in MultiQC config."
        )
    if "variants_file" not in plugin_config:
        raise ValueError(
            "Missing required 'multiqc_dumpling.variants_file' setting in MultiQC config."
        )

    # Parse the config file to get the orf length from the string
    # The orf is given in nucleotide coordinates in the form "XX-YY",
    # where XX and YY are integers.

    coords = plugin_config["orf"].split("-")
    # Check that the coordinates are specified in the expected format
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
    config.orf_length = end - start + 1

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
