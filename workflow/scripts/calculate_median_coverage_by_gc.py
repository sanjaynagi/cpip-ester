#!/usr/bin/env python3
"""
Calculate median coverage by GC content for CNV calling normalization
This script processes multiple samples at once, maintaining the 
functionality of the original script while adapting to the new input format
and Snakemake workflow.
"""
import sys
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np
import time
import socket
from re import sub

# Get inputs and outputs from Snakemake
gc_content_file = snakemake.input.gc_content
accessibility_files = snakemake.input.accessibility
mapq0_files = snakemake.input.mapq0
counts_files = snakemake.input.counts
accessibility_threshold = snakemake.params.accessibility_threshold
mapq0_threshold = snakemake.params.mapq0_threshold
median_coverage_output = snakemake.output.median_coverage
variance_output = snakemake.output.variance

metadata = pd.read_csv("results/metadata.tsv", sep="\t")
sample_names = metadata['Sample Name'].to_numpy()

print('Running script at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + 
      ' using machine ' + socket.gethostname() + '\n\n')
print('Input parameters:')
print(f'\taccessibility_threshold: {accessibility_threshold}')
print(f'\tgc_content_file: {gc_content_file}')
print(f'\tmapq0_threshold: {mapq0_threshold}')

# Function to normalize coverage by GC content (from original script)
def normalise_coverage_by_GC(coverage, median_coverage_by_GC, ploidy=2):
    """Normalize coverage values by GC content"""
    output = coverage.copy()
    # For each counts value in the output object, associate it with the median coverage for its GC bin
    output['expcov'] = [median_coverage_by_GC.loc[x] for x in output['GC']]
    # Now divide each counts value by the coverage associated with its GC bin
    output['Normalised_coverage'] = ploidy * output['Counts_total'] / output['expcov']
    # In rare cases, the expected coverage ends up being 0
    which_zeros = (output['expcov'] == 0)
    output.loc[(which_zeros, 'Normalised_coverage')] = 0
    return output

# Load GC content data
print("Loading GC content data...")
gc_data = pd.read_csv(gc_content_file, sep='\t')
# Rename columns to match expected format
gc_data = gc_data.rename(columns={
    'chrom': 'Chrom',
    'start': 'Position',
    'gc_content': 'GC'
})

# Take the floor of the percentage GC content
gc_data['GC'] = (gc_data['GC'] * 100).astype(int)

# Function to filter by quantile
def filter_by_quantile(x, thresh):
    return (x[x < np.quantile(x, thresh)])

# Initialize output tables
output_table = pd.DataFrame()
output_variance = pd.DataFrame(0, columns=[], index=sample_names)
chromosomes_all = []

# Process each sample
for i, sample in enumerate(sample_names):
    print(f"\nProcessing sample {sample} ({i+1}/{len(sample_names)})")
    
    # Load counts data for this sample
    counts_file = counts_files[i]
    print(f"Loading counts data from {counts_file}...")
    counts_data = pd.read_csv(counts_file, sep='\t')
    
    # Rename columns to match original script if needed
    if 'Chromosome' in counts_data.columns:
        counts_data = counts_data.rename(columns={'Chromosome':'Chrom'})
    if 'Counts total' in counts_data.columns:
        counts_data = counts_data.rename(columns={'Counts total':'Counts_total'})
    
    # Load accessibility data for this sample
    accessibility_file = accessibility_files[i]
    print(f"Loading accessibility data from {accessibility_file}...")
    accessibility_data = pd.read_csv(accessibility_file, sep='\t')
    # Apply accessibility threshold
    accessibility_data['filterpass'] = accessibility_data['Mean_accessibility'] >= accessibility_threshold
    
    # Load mapq0 data for this sample
    mapq0_file = mapq0_files[i]
    print(f"Loading mapq0 data from {mapq0_file}...")
    mapq0_data = pd.read_csv(mapq0_file, sep='\t')
    mapq0_data['prop_mapq0'] = mapq0_data['Count mapq = 0'] / (mapq0_data['Count mapq > 0'] + mapq0_data['Count mapq = 0'])
    mapq0_data['filterpass'] = mapq0_data.prop_mapq0 <= mapq0_threshold
    
    # Ensure all Position columns are numeric
    counts_data['Position'] = pd.to_numeric(counts_data['Position'])
    gc_data['Position'] = pd.to_numeric(gc_data['Position'])
    accessibility_data['Position'] = pd.to_numeric(accessibility_data['Position'])
    mapq0_data['Position'] = pd.to_numeric(mapq0_data['Position'])
    
    # Merge data
    print("Merging datasets...")
    data = pd.merge(counts_data, gc_data[['Chrom', 'Position', 'GC']], on=['Chrom', 'Position'])
    data = pd.merge(data, accessibility_data[['Chrom', 'Position', 'filterpass']], on=['Chrom', 'Position'])
    data = pd.merge(data, mapq0_data[['Chrom', 'Position', 'filterpass']], 
                    on=['Chrom', 'Position'], suffixes=('_acc', '_mapq0'))
    
    # Create combined filter pass
    data['doublefilterpass'] = data['filterpass_acc'] & data['filterpass_mapq0']
    # Only use data from accessible regions with low mapq0 proportion
    accessible_data = data[data['doublefilterpass']]
    
    # Identify autosomes and sex chromosomes
    chromosomes = accessible_data['Chrom'].unique()
    chromosomes_all.extend(chromosomes)
    sex_chroms = [c for c in chromosomes if c == 'X' or c == 'chrX']
    if len(sex_chroms) > 0:
        sex_chrom = sex_chroms[0]
    else:
        sex_chrom = None
    autosomes = [c for c in chromosomes if c not in sex_chroms]
    
    print(f"Identified autosomes: {autosomes}")
    if sex_chrom:
        print(f"Identified sex chromosome: {sex_chrom}")
    
    # Filter autosomal data
    autosomal_data = accessible_data[accessible_data['Chrom'].isin(autosomes)]
    
    # Group by GC content and calculate median coverage
    gc_bins = autosomal_data.groupby('GC')
    median_counts_by_GC = gc_bins['Counts_total'].median()
    
    # Create or update output dataframe for median coverage
    if i == 0:
        bin_counts = gc_bins.size().rename('bin_freq')
        output_table = pd.DataFrame({'bin_freq': bin_counts})
        output_table = output_table.loc[output_table.bin_freq > 0, :].copy()
        output_table.index.name = 'GC'
    
    # Add this sample's median coverage to the output table
    output_table[sample] = median_counts_by_GC
    
    # Normalize coverage for autosomes and sex chromosome
    print("Calculating normalized coverage...")
    autosomal_masked_counts = normalise_coverage_by_GC(autosomal_data, median_counts_by_GC)
    
    if sex_chrom:
        sexchrom_data = accessible_data[accessible_data['Chrom'] == sex_chrom]
        sexchrom_masked_counts = normalise_coverage_by_GC(sexchrom_data, median_counts_by_GC)
    
    # Calculate variance with filtering by quantile
    quantile_threshold = 0.99
    print("Calculating variance by chromosome...")
    
    # Calculate variance for autosomes
    autosomal_counts_by_chrom = autosomal_masked_counts.groupby('Chrom')
    for chrom in autosomes:
        if chrom not in output_variance.columns:
            output_variance[chrom] = np.nan
        chrom_data = autosomal_counts_by_chrom.get_group(chrom)
        output_variance.loc[sample, chrom] = np.var(
            filter_by_quantile(chrom_data.Normalised_coverage, quantile_threshold))
    
    # Calculate variance for sex chromosome if present
    if sex_chrom:
        if sex_chrom not in output_variance.columns:
            output_variance[sex_chrom] = np.nan
        output_variance.loc[sample, sex_chrom] = np.var(
            filter_by_quantile(sexchrom_masked_counts.Normalised_coverage, quantile_threshold))
    
    # Calculate variance across all autosomes
    if 'autosomes' not in output_variance.columns:
        output_variance['autosomes'] = np.nan
    output_variance.loc[sample, 'autosomes'] = np.var(
        filter_by_quantile(autosomal_masked_counts.Normalised_coverage, quantile_threshold))

# Ensure all chromosomes are in the variance dataframe
chromosomes_all = list(set(chromosomes_all))
for chrom in chromosomes_all:
    if chrom not in output_variance.columns:
        output_variance[chrom] = np.nan

# Save outputs with proper formatting
print(f"\nSaving median coverage output to file {median_coverage_output}")
output_table.to_csv(median_coverage_output, sep='\t')

print(f"\nSaving variance output to file {variance_output}")
output_variance.to_csv(variance_output, sep='\t')

print('\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
