#!/usr/bin/env python3
"""
Apply HMM to normalized coverage data for CNV calling
"""

import pandas as pd
import numpy as np
from hmmlearn.hmm import GaussianHMM

# Get inputs, outputs, and parameters from Snakemake
counts_file = snakemake.input.counts
gc_content_file = snakemake.input.gc_content
median_coverage_file = snakemake.input.median_coverage
variance_file = snakemake.input.variance
accessibility_file = snakemake.input.accessibility
mapq0_proportions_file = snakemake.input.mapq0_proportions
hmm_output_file = snakemake.output.hmm_output
cnv_calls_file = snakemake.output.cnv_calls

transition_probability = float(snakemake.params.transition_probability)
max_copy_number = int(snakemake.params.max_copy_number)
max_mapq0_proportion = float(snakemake.params.max_mapq0_proportion)

print(f"Running HMM CNV detection with transition probability: {transition_probability}")
print(f"Maximum copy number: {max_copy_number}")

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
    return output

# Load all required data
print("Loading input data...")
counts_data = pd.read_csv(counts_file, sep='\t')
gc_data = pd.read_csv(gc_content_file, sep='\t', names=['Chrom', 'Position', 'GC'])
median_coverage_by_gc = pd.read_csv(median_coverage_file, sep='\t', index_col='GC')
variance_data = pd.read_csv(variance_file, sep='\t')
accessibility_data = pd.read_csv(accessibility_file, sep='\t')
mapq0_data = pd.read_csv(mapq0_proportions_file, sep='\t')

# Create MAPQ0 mask
mapq0_data['mapq0_acceptable'] = mapq0_data['prop_mapq0'] <= max_mapq0_proportion

# Merge data for processing
data = pd.merge(counts_data, gc_data, on=['Chrom', 'Position'])
data = pd.merge(data, mapq0_data[['Chrom', 'Position', 'mapq0_acceptable']], on=['Chrom', 'Position'])
data = pd.merge(data, accessibility_data[['Chrom', 'Position', 'filterpass']], on=['Chrom', 'Position'])

# Apply masks - both MAPQ0 and accessibility
print("Applying filters...")
masked_data = data[(data['mapq0_acceptable']) & (data['filterpass'])]

# Only use GC bins with sufficient data
valid_gc_bins = median_coverage_by_gc[median_coverage_by_gc['bin_freq'] >= 100].index
masked_data = masked_data[masked_data['GC'].isin(valid_gc_bins)]

print(f"After filtering, {len(masked_data)} windows remain for analysis")

# Normalize coverage by GC content
print("Normalizing coverage...")
normalized_data = normalize_coverage_by_gc(masked_data, median_coverage_by_gc['median_coverage'])

# Get autosomal variance
autosomal_variance = variance_data['autosomes'].iloc[0]

# Run HMM
print("Running HMM for CNV detection...")
cnv_states = fit_hmm(
    normalized_data['Normalized_coverage'].values,
    transition_probability=transition_probability,
    variance=autosomal_variance/2,  # We divide by 2 as the variance is for diploid state
    max_copy_number=max_copy_number
)

# Add CNV states to the data
normalized_data['CNV'] = cnv_states

# Save HMM output
normalized_data.to_csv(hmm_output_file, sep='\t', index=False)

# Process CNV calls - merge adjacent windows with the same copy number
print("Processing CNV calls...")
normalized_data['CNV_change'] = normalized_data['CNV'].diff() != 0
cnv_segments = []

current_segment = None
for _, row in normalized_data.iterrows():
    if row['CNV_change'] or current_segment is None:
        if current_segment is not None:
            cnv_segments.append(current_segment)
        
        current_segment = {
            'Chromosome': row['Chrom'],
            'Start': row['Position'],
            'End': row['Position'] + int(snakemake.params.window_size),
            'CopyNumber': row['CNV'],
            'MeanNormalizedCoverage': row['Normalized_coverage']
        }
    else:
        current_segment['End'] = row['Position'] + int(snakemake.params.window_size)
        current_segment['MeanNormalizedCoverage'] = (current_segment['MeanNormalizedCoverage'] + 
                                                    row['Normalized_coverage']) / 2

# Add the last segment
if current_segment is not None:
    cnv_segments.append(current_segment)

# Create and save CNV calls dataframe
cnv_calls = pd.DataFrame(cnv_segments)
cnv_calls.to_csv(cnv_calls_file, index=False)

print(f"HMM processing complete. Outputs saved to {hmm_output_file} and {cnv_calls_file}")