#!/usr/bin/env python3
"""
Calculate the proportion of MAPQ0 reads in genomic windows
"""
import sys
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np
import pysam

# Get inputs, outputs, and parameters from Snakemake
bam_file = snakemake.input.bam
output_file = snakemake.output.mapq0_proportions
window_size = int(snakemake.params.window_size)
interval = int(snakemake.params.interval)

print(f"Calculating MAPQ0 proportions for: {bam_file}")
print(f"Window size: {window_size}bp, Interval: {interval}bp")

# Open the BAM file
alignment = pysam.AlignmentFile(bam_file, 'rb')

# Get chromosomes and their lengths
chromosomes = alignment.references
chromosome_lengths = alignment.lengths

# Initialize results dataframe
results = []

# Process each chromosome
for chrom, chrom_length in zip(chromosomes, chromosome_lengths):
    print(f"Processing chromosome {chrom} ({chrom_length}bp)")
    
    # Define windows
    positions = range(0, chrom_length, interval)
    
    # Process each window
    for pos in positions:
        # Fetch alignments in this window
        alignments = alignment.fetch(chrom, pos, min(pos + window_size, chrom_length))
        
        # Count reads by MAPQ
        count_mapq0 = 0
        count_mapq_above0 = 0
        
        for read in alignments:
            if not read.is_unmapped:
                if read.reference_start >= pos:
                    if read.mapping_quality == 0:
                        count_mapq0 += 1
                    else:
                        count_mapq_above0 += 1
        
        # Calculate proportion
        total_reads = count_mapq0 + count_mapq_above0
        prop_mapq0 = count_mapq0 / total_reads if total_reads > 0 else np.nan
        
        # Store results
        results.append({
            'Chrom': chrom,
            'Position': pos,
            'Count mapq = 0': count_mapq0,
            'Count mapq > 0': count_mapq_above0,
            'prop_mapq0': prop_mapq0
        })

# Create dataframe and save results
df = pd.DataFrame(results)
df.to_csv(output_file, sep='\t', index=False)

print(f"Completed processing. Output saved to {output_file}")