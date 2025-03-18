#!/usr/bin/env python3
"""
Calculate genomic accessibility based on GC content and MAPQ0 proportions
"""

import pandas as pd
import numpy as np

# Get inputs, outputs, and parameters from Snakemake
gc_content_file = snakemake.input.gc_content
mapq0_proportions_file = snakemake.output.mapq0_proportions
output_file = snakemake.output.accessibility
mapq0_threshold = float(snakemake.params.mapq0_threshold)
accessibility_threshold = float(snakemake.params.accessibility_threshold)

print(f"Calculating accessibility using GC content: {gc_content_file}")
print(f"MAPQ0 threshold: {mapq0_threshold}, Accessibility threshold: {accessibility_threshold}")

# Load GC content data
gc_data = pd.read_csv(gc_content_file, sep='\t', names=['Chrom', 'Position', 'GC'])

# Load MAPQ0 proportions data
mapq0_data = pd.read_csv(mapq0_proportions_file, sep='\t')

# Make sure both dataframes have the same number of rows
if len(gc_data) != len(mapq0_data):
    raise ValueError("GC content and MAPQ0 proportion files have different number of rows")

# Create accessibility mask based on MAPQ0 proportions
mapq0_data['filterpass'] = mapq0_data['prop_mapq0'] <= mapq0_threshold

# Calculate accessibility (this would typically involve more complex criteria,
# but for simplicity we're using MAPQ0 proportions as the main criterion)
results = pd.DataFrame({
    'Chrom': gc_data['Chrom'],
    'Position': gc_data['Position'],
    'GC': gc_data['GC'],
    'Mean_accessibility': np.where(mapq0_data['filterpass'], 1.0, 0.0),
    'filterpass': mapq0_data['filterpass']
})

# Save to file
results.to_csv(output_file, sep='\t', index=False)

print(f"Completed accessibility calculation. Output saved to {output_file}")