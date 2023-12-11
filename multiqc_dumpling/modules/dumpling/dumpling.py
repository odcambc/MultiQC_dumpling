#!/usr/bin/env python

from __future__ import print_function

import logging
from collections import OrderedDict

import numpy as np
import pandas as pd

from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")
log.info("Loaded MultiQC_dumpling v%s", get_distribution("multiqc_dumpling").version)


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

        # Filter out samples matching ignored sample names
        self.dumpling_count_plot_data = self.ignore_samples(
            self.dumpling_count_plot_data
        )
        self.dumpling_count_data = self.ignore_samples(self.dumpling_count_data)

        self.dumpling_coverage_data = self.ignore_samples(self.dumpling_coverage_data)
        self.dumpling_coverage_plot_data = self.ignore_samples(
            self.dumpling_coverage_plot_data
        )

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.dumpling_count_plot_data) == 0:
            log.debug("Could not find any reports in %s", config.analysis_dir)
            raise ModuleNoSamplesFound

        log.info("Found %s processed count files", len(self.dumpling_count_plot_data))

        if len(self.dumpling_coverage_data) == 0:
            log.debug("Could not find any coverage reports in %s", config.analysis_dir)
            raise ModuleNoSamplesFound

        log.info("Found %s coverage report files", len(self.dumpling_coverage_data))

        # Write all the parsed report data to files
        self.write_data_file(self.dumpling_count_data, "multiqc_dumpling_counts")
        self.write_data_file(self.dumpling_coverage_data, "multiqc_dumpling_coverage")
        self.write_data_file(self.dumpling_count_plot_data, "multiqc_dumpling_counts_plot")
        self.write_data_file(self.dumpling_coverage_plot_data, "multiqc_dumpling_coverage_plot")

        # Add to the table
        headers = OrderedDict()
        headers = {
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
            },
        }

        self.general_stats_addcols(self.dumpling_count_data, headers)
        self.general_stats_addcols(self.dumpling_coverage_data, headers)

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

    def parse_counts(self, file, bin_n=30):
        """Parse the counts file and return a dict of counts and their frequencies for
        plotting a histogram.

        At the same time, set the variants number if it doesn't already exist.
        """
        df = pd.read_csv(file, index_col=False, header=0)

        max_counts = int(df["count"].max())
        mean_counts = df["count"].mean()
        median_counts = df["count"].median()
        n_zero_counts = int(len(df[df["count"] == 0]))

        if not hasattr(config, "n_variants"):
            log.warning(
                "Number of variants not set: using counts file %s, setting n_variants to %s",
                file,
                len(df),
            )
            config.n_variants = len(df)

        fraction_zero_counts = n_zero_counts / config.n_variants

        # If there are no counts, return a dict with a single entry.

        if max_counts == 0:
            return {0: len(df)}

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

        counts_stats_dict = {
            "Max counts": max_counts,
            "Mean counts": mean_counts,
            "Median counts": median_counts,
            "Number of zero counts": n_zero_counts,
            "Fraction of zero counts": fraction_zero_counts,
        }
        
        return counts_bin_dict, counts_stats_dict
        
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

        max_coverage = df["Coverage"].max()
        min_coverage = df["Coverage"].min()
        mean_coverage = df["Coverage"].mean()

        if max_coverage < bin_n:
            bin_n = max_coverage


        bin_n = max(1, int(max_coverage))

        logging.debug("Using %x bins", bin_n)

        out, bins = pd.cut(
            df["Coverage"], bins=bin_n, include_lowest=True, right=False, retbins=True
        )

        # now make dict between bin and out

        coverage_bin_dict = dict(
            zip(bins.tolist(), out.value_counts(normalize=False).sort_index().tolist())
        )

        #coverage_bin_dict = df["Coverage"].value_counts(normalize=False).to_dict()

        coverage_stats_dict = {
            "Max coverage": max_coverage,
            "Min coverage": min_coverage,
            "Mean coverage": mean_coverage,
        }

        return coverage_bin_dict, coverage_stats_dict