#!/usr/bin/env python3
"""
Apply HMM to normalized coverage data for CNV calling
Process multiple samples in a single run
"""
import sys
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np
from hmmlearn.hmm import GaussianHMM
import os
import time
import socket
from sys import stdout

# Get inputs, outputs, and parameters from Snakemake
counts_files_dir = snakemake.params.counts_dir  # Directory with all count files
gc_content_file = snakemake.input.gc_content
median_coverage_file = snakemake.input.median_coverage  # Now contains data for all samples
variance_file = snakemake.input.variance  # Now contains data for all samples
accessibility_file = snakemake.input.accessibility
mapq0_proportions_file = snakemake.input.mapq0_proportions
hmm_output_dir = snakemake.params.hmm_output_dir  # Directory for output files

metadata = pd.read_csv("results/metadata.tsv", sep="\t")
sample_names = metadata['Sample Name'].to_numpy()

transition_probability = float(snakemake.params.transition_probability)
max_copy_number = int(snakemake.params.max_copy_number)
max_mapq0_proportion = float(snakemake.params.max_mapq0_proportion)
window_size = int(snakemake.params.window_size)

print('Running script at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + 
      ' using machine ' + socket.gethostname() + '\n\n')
print(f"Running HMM CNV detection with transition probability: {transition_probability}")
print(f"Maximum copy number: {max_copy_number}")
print(f"Maximum MAPQ0 proportion: {max_mapq0_proportion}")
print(f"Processing {len(sample_names)} samples")


def fit_hmm(depth_normed, transition_probability, variance, variance_fixed=0.001, 
           max_copy_number=12, n_iter=0, params='st', init_params=''):
    """
    Fit an HMM to normalized depth data to infer copy number states
    
    Parameters:
    depth_normed: Normalized depth values
    transition_probability: Probability of state transitions
    variance: Variance per copy number
    variance_fixed: Fixed variance component for zero copy state
    max_copy_number: Maximum copy number to consider
    n_iter: Number of iterations for HMM fitting
    params: Parameters to estimate during fitting
    init_params: Parameters to initialize from data
    
    Returns:
    Array of copy number states for each position
    """
    # Set up copy number range
    min_copy_number = 0
    n_states = max_copy_number - min_copy_number + 1
    
    # Construct the transition matrix
    transmat = np.zeros((n_states, n_states))
    transmat[:] = transition_probability
    transmat[np.diag_indices(n_states)] = 1 - ((n_states - 1) * transition_probability)
    
    # Construct means and covariance
    means_list = list(range(n_states))
    means = np.array([[n] for n in means_list])
    covars = np.array([[variance * n + variance_fixed] for n in means_list])
    
    # Setup HMM
    model = GaussianHMM(n_states,
                      covariance_type='diag',
                      n_iter=n_iter,
                      params=params,
                      init_params=init_params)
    model.means_ = means
    model.covars_ = covars
    model.transmat_ = transmat
    
    # Fit HMM and predict states
    obs = np.column_stack([depth_normed])
    model.fit(obs)
    hidden_states = model.predict(obs)
    
    return hidden_states

def normalize_coverage_by_gc(coverage, median_coverage_by_gc, ploidy=2):
    """
    Normalize coverage data by GC content
    
    Parameters:
    coverage: DataFrame with coverage data
    median_coverage_by_gc: Series with median coverage by GC content
    ploidy: Expected genome ploidy (typically 2)
    
    Returns:
    DataFrame with normalized coverage values
    """
    output = coverage.copy()
    output['Expected_coverage'] = [median_coverage_by_gc.loc[x] if x in median_coverage_by_gc.index else np.nan 
                                  for x in output['GC']]
    output['Normalized_coverage'] = ploidy * output['Counts_total'] / output['Expected_coverage']
    
    # Handle cases where expected coverage is 0 or NaN
    output.loc[output['Expected_coverage'] == 0, 'Normalized_coverage'] = 0
    output.loc[np.isnan(output['Expected_coverage']), 'Normalized_coverage'] = np.nan
    
    return output

# Load shared reference data
print("Loading reference data...")


# Load GC data with correct column names
gc_data = pd.read_csv(gc_content_file, sep='\t')
# Rename columns to match expected format
gc_data = gc_data.rename(columns={
    'chrom': 'Chrom',
    'start': 'Position',
    'gc_content': 'GC'
})
# Take the floor of the percentage GC content
gc_data['GC'] = (gc_data['GC'] * 100).astype(int)

# Load median coverage for all samples
median_coverage_by_gc = pd.read_csv(median_coverage_file, sep='\t', index_col='GC')

# Load variance data for all samples
variance_data = pd.read_csv(variance_file, sep='\t', index_col=0)

# Ensure numeric position columns
gc_data['Position'] = pd.to_numeric(gc_data['Position'])

# Get valid GC bins with sufficient data
valid_gc_bins = median_coverage_by_gc[median_coverage_by_gc['bin_freq'] >= 100].index

# Process each sample
for sample_idx, sample in enumerate(sample_names):
    print(f"\nProcessing sample {sample} ({sample_idx+1}/{len(sample_names)})...")
    
    # Define output files for this sample
    hmm_output_file = os.path.join(hmm_output_dir, sample, "hmm_output.tsv")
    cnv_calls_file = os.path.join(hmm_output_dir, sample, "cnv_calls.csv")
    
    # Load this sample's count data
    counts_file = os.path.join(counts_files_dir, sample, "counts_for_hmm.tsv")
    counts_data = pd.read_csv(counts_file, sep='\t')

    # Load accessibililty and MAPQ0 data
    accessibility_data = pd.read_csv(f"results/cnv/{sample}/accessibility.tsv", sep='\t')
    mapq0_data = pd.read_csv(f"results/cnv/{sample}/mapq0_proportions.tsv", sep='\t')
    accessibility_data['Position'] = pd.to_numeric(accessibility_data['Position'])
    mapq0_data['Position'] = pd.to_numeric(mapq0_data['Position'])
    mapq0_data['mapq0_acceptable'] = mapq0_data['prop_mapq0'] <= max_mapq0_proportion
    
    # Rename columns if needed
    if 'Chromosome' in counts_data.columns:
        counts_data = counts_data.rename(columns={'Chromosome':'Chrom'})
    if 'Counts total' in counts_data.columns:
        counts_data = counts_data.rename(columns={'Counts total':'Counts_total'})
    
    counts_data['Position'] = pd.to_numeric(counts_data['Position'])
    
    # Merge data for processing
    print(f"  Merging datasets for sample {sample}...")
    data = pd.merge(counts_data, gc_data[['Chrom', 'Position', 'GC']], on=['Chrom', 'Position'])
    data = pd.merge(data, mapq0_data[['Chrom', 'Position', 'mapq0_acceptable']], on=['Chrom', 'Position'])
    data = pd.merge(data, accessibility_data[['Chrom', 'Position', 'filterpass']], on=['Chrom', 'Position'])
    
    # Apply masks - both MAPQ0 and accessibility
    print(f"  Applying filters for sample {sample}...")
    masked_data = data[(data['mapq0_acceptable']) & (data['filterpass'])]
    
    # Only use GC bins with sufficient data
    masked_data = masked_data[masked_data['GC'].isin(valid_gc_bins)]
    
    print(f"  After filtering, {len(masked_data)} windows remain for analysis of sample {sample}")
    if len(masked_data) == 0:
        print(f"Warning: No data remains after filtering for sample {sample}. Skipping HMM.")
        continue
    
    # Get the median coverage values specific to this sample
    if sample in median_coverage_by_gc.columns:
        sample_median_coverage = median_coverage_by_gc[sample]
    else:
        print(f"Warning: No median coverage data for sample {sample}. Using first available sample.")
        # Use first sample's data as fallback
        first_sample_col = next(col for col in median_coverage_by_gc.columns if col != 'bin_freq')
        sample_median_coverage = median_coverage_by_gc[first_sample_col]
    
    # Normalize coverage by GC content
    print(f"  Normalizing coverage for sample {sample}...")
    normalized_data = normalize_coverage_by_gc(masked_data, sample_median_coverage)
    
    # Filter out NaN values
    normalized_data = normalized_data.dropna(subset=['Normalized_coverage'])
    
    # Get this sample's autosomal variance
    if sample in variance_data.index and 'autosomes' in variance_data.columns:
        autosomal_variance = variance_data.loc[sample, 'autosomes']
    else:
        print(f"Warning: No variance data for sample {sample}. Using average variance.")
        # Use average of available variances
        autosomal_variance = variance_data['autosomes'].mean()
    
    print(f"  Sample {sample} autosomal variance: {autosomal_variance}")
    print(f"  Running HMM for CNV detection on sample {sample}...")
    cnv_states = fit_hmm(
        normalized_data['Normalized_coverage'].values,
        transition_probability=transition_probability,
        variance=autosomal_variance/2,  # We divide by 2 as the variance is for diploid state
        max_copy_number=max_copy_number
    )
    
    # Add CNV states to the data
    normalized_data['CNV'] = cnv_states
    
    # Save HMM output
    print(f"  Saving HMM output for sample {sample}...")
    normalized_data.to_csv(hmm_output_file, sep='\t', index=False)
    
    # Process CNV calls - merge adjacent windows with the same copy number
    print(f"  Processing CNV calls for sample {sample}...")
    normalized_data = normalized_data.sort_values(['Chrom', 'Position'])
    normalized_data['CNV_change'] = (normalized_data['CNV'].diff() != 0) | (normalized_data['Chrom'].shift() != normalized_data['Chrom'])
    cnv_segments = []
    
    current_segment = None
    for _, row in normalized_data.iterrows():
        if row['CNV_change'] or current_segment is None:
            if current_segment is not None:
                cnv_segments.append(current_segment)
            
            current_segment = {
                'Sample': sample,
                'Chromosome': row['Chrom'],
                'Start': row['Position'],
                'End': row['Position'] + window_size,
                'CopyNumber': row['CNV'],
                'MeanNormalizedCoverage': row['Normalized_coverage'],
                'WindowCount': 1
            }
        else:
            current_segment['End'] = row['Position'] + window_size
            # Running average
            current_segment['MeanNormalizedCoverage'] = (
                (current_segment['MeanNormalizedCoverage'] * current_segment['WindowCount']) + 
                row['Normalized_coverage']
            ) / (current_segment['WindowCount'] + 1)
            current_segment['WindowCount'] += 1
    
    # Add the last segment
    if current_segment is not None:
        cnv_segments.append(current_segment)
    
    # Create and save CNV calls dataframe
    cnv_calls = pd.DataFrame(cnv_segments)
    cnv_calls.to_csv(cnv_calls_file, index=False)
    
    print(f"  Completed processing for sample {sample}")

print("\nHMM processing complete for all samples.")
print('Script finished at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')