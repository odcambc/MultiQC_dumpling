#!/usr/bin/env python

from __future__ import print_function
import importlib_metadata

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

    # Parse the config file to get the orf length from the string
    # The orf is given in nucleotide coordinates in the form "XX-YY",
    # where XX and YY are integers.

    coords = config.multiqc_dumpling['orf'].split('-')
    # Check that the coordinates are specified in the expected format
    if len(coords) != 2:
        log.error("Invalid ORF coordinates: %s", config.multiqc_dumpling['orf'])
        raise ValueError("Invalid ORF coordinates")
    if not all(coord.isdigit() for coord in coords):
        log.error("Invalid ORF coordinates: %s", config.multiqc_dumpling['orf'])
        raise ValueError("Invalid ORF coordinates")
    if int(coords[0]) > int(coords[1]):
        log.error("Invalid ORF coordinates: %s", config.multiqc_dumpling['orf'])
        raise ValueError("Invalid ORF coordinates")
    if int(coords[0]) < 1:
        log.error("Invalid ORF coordinates: %s", config.multiqc_dumpling['orf'])
        raise ValueError("Invalid ORF coordinates")
    start = int(coords[0])
    end = int(coords[1])
    config.orf_length = end - start + 1

    # Read the designed variants file to get the number of designed variants
    # The first line is a header, so subtract 1 from the total number of lines

    with open(config.multiqc_dumpling['variants_file'], "r") as f:
        config.n_variants = sum(1 for line in f) - 1

    # Add the designed variants file to the ignore list
    config.fn_ignore_files.extend([config.multiqc_dumpling['variants_file']])