#!/usr/bin/env python3
"""
Calculate median coverage by GC content for CNV calling normalization
"""

import pandas as pd
import numpy as np

# Get inputs and outputs from Snakemake
counts_file = snakemake.input.counts
gc_content_file = snakemake.input.gc_content
accessibility_file = snakemake.input.accessibility
median_coverage_output = snakemake.output.median_coverage
variance_output = snakemake.output.variance

print(f"Calculating median coverage by GC using counts from: {counts_file}")

# Load input data
counts_data = pd.read_csv(counts_file, sep='\t')
gc_data = pd.read_csv(gc_content_file, sep='\t', names=['Chrom', 'Position', 'GC'])
accessibility_data = pd.read_csv(accessibility_file, sep='\t')

# Merge data
data = pd.merge(counts_data, gc_data, on=['Chrom', 'Position'])
data = pd.merge(data, accessibility_data[['Chrom', 'Position', 'filterpass']], on=['Chrom', 'Position'])

# Only use data from accessible regions
accessible_data = data[data['filterpass']]

# Group by GC content and calculate median coverage
gc_bins = accessible_data.groupby('GC')
median_coverage_by_gc = gc_bins['Counts_total'].median()

# Create output dataframe for median coverage
bin_counts = gc_bins.size().rename('bin_freq')
median_coverage_output_df = pd.DataFrame({'median_coverage': median_coverage_by_gc, 'bin_freq': bin_counts})
median_coverage_output_df.index.name = 'GC'

# Save median coverage by GC content
median_coverage_output_df.to_csv(median_coverage_output, sep='\t')

# Function to normalize coverage by GC content
def normalize_coverage_by_gc(coverage_data, median_by_gc, ploidy=2):
    """Normalize coverage values by GC content"""
    normalized = coverage_data.copy()
    normalized['expected_coverage'] = normalized['GC'].map(median_by_gc)
    normalized['normalized_coverage'] = ploidy * normalized['Counts_total'] / normalized['expected_coverage']
    
    # Handle cases where expected coverage is 0
    normalized.loc[normalized['expected_coverage'] == 0, 'normalized_coverage'] = 0
    
    return normalized

# Normalize coverage
normalized_data = normalize_coverage_by_gc(accessible_data, median_coverage_by_gc)

# Calculate variance
def calculate_variance(data, quantile_threshold=0.99):
    """Calculate variance after filtering extreme values"""
    filtered_data = data[data < np.quantile(data, quantile_threshold)]
    return np.var(filtered_data)

# Calculate variance by chromosome
chromosomes = accessible_data['Chrom'].unique()
autosomes = [c for c in chromosomes if not (c.startswith('X') or c.startswith('Y') or c.startswith('chrX') or c.startswith('chrY'))]

# Create dictionary for variance values
variance_dict = {'autosomes': calculate_variance(normalized_data['normalized_coverage'])}

for chrom in chromosomes:
    chrom_data = normalized_data[normalized_data['Chrom'] == chrom]
    variance_dict[chrom] = calculate_variance(chrom_data['normalized_coverage'])

# Save variance data
variance_df = pd.DataFrame([variance_dict])
variance_df.to_csv(variance_output, sep='\t', index=False)

print(f"Completed calculation. Outputs saved to {median_coverage_output} and {variance_output}")