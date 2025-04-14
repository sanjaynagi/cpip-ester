#!/usr/bin/env python3
"""
Count reads in genomic windows for HMM-based CNV detection
"""
import sys
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np
import pysam

# Get inputs, outputs, and parameters from Snakemake
bam_file = snakemake.input.bam
output_file = snakemake.output.counts
window_size = int(snakemake.params.window_size)
interval = int(snakemake.params.interval)
min_mapq = int(snakemake.params.min_mapq)

print(f"Processing BAM file: {bam_file}")
print(f"Window size: {window_size}bp, Interval: {interval}bp, Min MAPQ: {min_mapq}")

# Open the BAM file
alignment = pysam.AlignmentFile(bam_file, 'rb')

# Get chromosomes and their lengths
chromosomes = np.array(alignment.references)
chromosome_lengths = np.array(alignment.lengths)

mask = chromosomes == "CM027412.1"
chrom = chromosomes[mask][0]
chrom_length = chromosome_lengths[mask][0]

# Initialize results dataframe
results = []
print(f"Processing chromosome {chrom} ({chrom_length}bp)")

# Define windows
positions = range(0, chrom_length, interval)
# Process each window
for pos in positions:
    # Fetch alignments in this window
    alignments = alignment.fetch(chrom, pos, min(pos + window_size, chrom_length))
    
    # Count reads with sufficient mapping quality
    count = 0
    count_mapq0 = 0
    count_below_threshold = 0
    
    for read in alignments:
        if not read.is_unmapped:
            # Count reads that start in the current window
            if read.reference_start >= pos:
                if read.mapping_quality >= min_mapq:
                    count += 1
                elif read.mapping_quality == 0:
                    count_mapq0 += 1
                else:
                    count_below_threshold += 1
    
    # Store results
    results.append({
        'Chrom': chrom,
        'Position': pos,
        'Counts_mapq_above_threshold': count,
        'Counts_below_threshold': count_below_threshold,
        'Counts_mapq0': count_mapq0,
        'Counts_total': count + count_mapq0 + count_below_threshold
    })

# Create dataframe and save results
df = pd.DataFrame(results)
df.to_csv(output_file, sep='\t', index=False)

print(f"Completed processing. Output saved to {output_file}")