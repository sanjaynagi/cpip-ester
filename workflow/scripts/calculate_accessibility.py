#!/usr/bin/env python3
"""
Calculate genomic accessibility based on GC content and MAPQ0 proportions
"""

import pandas as pd
import numpy as np

# Get inputs, outputs, and parameters from Snakemake
gc_content_file = snakemake.input.gc_content
mapq0_proportions_file = snakemake.input.mapq0_proportions
output_file = snakemake.output.accessibility
mapq0_threshold = float(snakemake.params.mapq0_threshold)
accessibility_threshold = float(snakemake.params.accessibility_threshold)

print(f"Calculating accessibility using GC content: {gc_content_file}")
print(f"MAPQ0 threshold: {mapq0_threshold}, Accessibility threshold: {accessibility_threshold}")

# Load GC content data with correct column names
gc_data = pd.read_csv(gc_content_file, sep='\t')

# Convert to expected format
gc_data = gc_data.rename(columns={
    'chrom': 'Chrom',
    'start': 'Position',  # Using start position as the Position
    'gc_content': 'GC'
})

# Load MAPQ0 proportions data
mapq0_data = pd.read_csv(mapq0_proportions_file, sep='\t')

# Make sure both dataframes have the same number of rows
if len(gc_data) != len(mapq0_data):
    raise ValueError("GC content and MAPQ0 proportion files have different number of rows")

# Create accessibility mask based on MAPQ0 proportions
mapq0_data['filterpass'] = mapq0_data['prop_mapq0'] <= mapq0_threshold

# First create the accessibility array
accessibility = np.where(mapq0_data['filterpass'], 1.0, 0.0)
gc_data = gc_data.reset_index(drop=True)
mapq0_data = mapq0_data.reset_index(drop=True)

# Then create the DataFrame with explicit index
results = pd.DataFrame({
    'Chrom': gc_data['Chrom'].values,
    'Position': gc_data['Position'].values,
    'GC': gc_data['GC'].values,
    'Mean_accessibility': accessibility,
    'filterpass': mapq0_data['filterpass'].values
})

# Include the end position in the results if needed
if 'end' in gc_data.columns:
    results['End'] = gc_data['end'].values

# Save to file
results.to_csv(output_file, sep='\t', index=False)

print(f"Completed accessibility calculation. Output saved to {output_file}")
