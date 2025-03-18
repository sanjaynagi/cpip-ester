#!/usr/bin/env python3
"""
Calculate modal CNVs for all genes in the GFF file
"""

import pandas as pd
import numpy as np
import gffutils
import tempfile
import os
from collections import Counter

# Get inputs, outputs, and parameters from Snakemake
cnv_files = snakemake.input.cnv_files
annotation_file = snakemake.input.annotation
modal_cnvs_file = snakemake.output.modal_cnvs
samples = snakemake.params.samples

# Create a GFF database for querying
print("Creating GFF database...")
db_file = tempfile.NamedTemporaryFile(delete=False).name
db = gffutils.create_db(
    annotation_file,
    dbfn=db_file,
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True,
    disable_infer_transcripts=True,
    disable_infer_genes=True
)

# Get all genes from the GFF file
print("Extracting all genes from GFF file...")
genes = list(db.features_of_type('gene'))
print(f"Found {len(genes)} genes in the GFF file")

def get_gene_name(gene):
    """
    Get the name of a gene feature
    
    Parameters:
    gene: Gene feature from gffutils
    
    Returns:
    String with gene name
    """
    # Try different attribute names for gene name
    for attr in ['ID', 'Name', 'gene_id', 'gene_name']:
        if attr in gene.attributes:
            return gene[attr][0]
    # If no name attributes found, use coordinates
    return f"{gene.seqid}:{gene.start}-{gene.end}"

def calculate_modal_copy_number(cnv_data, chrom, start, end):
    """
    Calculate the modal copy number for a region
    
    Parameters:
    cnv_data: DataFrame with CNV calls
    chrom: Chromosome name
    start: Start position
    end: End position
    
    Returns:
    The most common copy number in the region
    """
    # Filter CNV data for the region
    region_cnvs = cnv_data[
        (cnv_data['Chromosome'] == chrom) & 
        (cnv_data['End'] > start) & 
        (cnv_data['Start'] < end)
    ]
    
    if len(region_cnvs) == 0:
        return 2  # Assume diploid if no CNVs in the region
    
    # For overlapping CNVs, calculate weighted copy number by overlap size
    weighted_copy_numbers = []
    
    for _, cnv in region_cnvs.iterrows():
        # Calculate overlap
        overlap_start = max(start, cnv['Start'])
        overlap_end = min(end, cnv['End'])
        overlap_size = overlap_end - overlap_start
        
        # Only count if there's actual overlap
        if overlap_size > 0:
            weighted_copy_numbers.extend([cnv['CopyNumber']] * overlap_size)
    
    if not weighted_copy_numbers:
        return 2  # Assume diploid if no overlap
    
    # Get modal copy number
    copy_number_counts = Counter(weighted_copy_numbers)
    modal_cn = copy_number_counts.most_common(1)[0][0]
    
    return modal_cn

# Process each CNV file and gene
results = []

for i, cnv_file in enumerate(cnv_files):
    sample = samples[i]
    print(f"Processing sample: {sample}")
    
    # Load CNV data
    cnv_data = pd.read_csv(cnv_file)
    
    # Process each gene
    for gene in genes:
        gene_name = get_gene_name(gene)
        chrom = gene.seqid
        start = gene.start
        end = gene.end
        
        # Calculate modal copy number
        modal_cn = calculate_modal_copy_number(cnv_data, chrom, start, end)
        
        # Store result
        results.append({
            'Sample': sample,
            'Region': gene_name,
            'Chromosome': chrom,
            'Start': start,
            'End': end,
            'ModalCopyNumber': modal_cn,
            'Length': end - start
        })

# Create and save results
results_df = pd.DataFrame(results)
results_df.to_csv(modal_cnvs_file, index=False)

# Clean up temporary database file
os.unlink(db_file)

print(f"Modal CNV analysis complete. Results saved to {modal_cnvs_file}")