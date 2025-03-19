#!/usr/bin/env python3
"""
Calculate GC content in genomic windows.
This script reads a reference genome, creates windows of specified size and step,
and calculates the GC content for each window.

To be used with Snakemake's script directive.
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from pyfaidx import Fasta
import pandas as pd


def setup_logging(log_file):
    """Set up logging to file and console."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def get_chromosome_sizes(fasta_file):
    """Get chromosome sizes from FASTA file."""
    fasta = Fasta(fasta_file)
    sizes = {name: len(seq) for name, seq in fasta.items()}
    return sizes


def create_windows(chrom, size, window_size, interval):
    """Create windows for a chromosome."""
    starts = np.arange(0, size - window_size + 1, interval)
    ends = starts + window_size
    # Make sure the last window doesn't exceed chromosome length
    if ends[-1] > size:
        ends[-1] = size
    return list(zip([chrom] * len(starts), starts, ends))


def calculate_gc_content(args):
    """Calculate GC content for a single window."""
    chrom, start, end = args
    try:
        sequence = fasta[chrom][start:end].seq.upper()
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        gc_count = g_count + c_count
        gc_content = gc_count / len(sequence) if len(sequence) > 0 else 0
        return chrom, start, end, gc_content
    except Exception as e:
        logging.error(f"Error processing {chrom}:{start}-{end}: {str(e)}")
        return chrom, start, end, None


# Get Snakemake parameters
reference_file = snakemake.input.reference
output_file = snakemake.output.gc_content
log_file = snakemake.log[0]
window_size = snakemake.params.window_size
interval = snakemake.params.interval
threads = snakemake.threads

# Set up logging
setup_logging(log_file)
logging.info(f"Starting GC content calculation with window size {window_size} and step {interval}")

# Load reference genome
global fasta
fasta = Fasta(reference_file)

# Get chromosome sizes
logging.info("Getting chromosome sizes")
chrom_sizes = get_chromosome_sizes(reference_file)

# Create windows for all chromosomes
logging.info("Creating windows")
all_windows = []
for chrom, size in chrom_sizes.items():
    windows = create_windows(chrom, size, window_size, interval)
    all_windows.extend(windows)

logging.info(f"Created {len(all_windows)} windows across {len(chrom_sizes)} chromosomes")

# Calculate GC content in parallel
logging.info(f"Calculating GC content using {threads} threads")
results = []

with ProcessPoolExecutor(max_workers=threads) as executor:
    for result in executor.map(calculate_gc_content, all_windows):
        if result[3] is not None:  # Skip failed windows
            results.append(result)

# Create dataframe and save to file
logging.info("Saving results to file")
df = pd.DataFrame(results, columns=["chrom", "start", "end", "gc_content"])
df.to_csv(output_file, sep="\t", index=False)

logging.info("Finished GC content calculation")