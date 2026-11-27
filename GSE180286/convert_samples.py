#!/usr/bin/env python3
"""
GSE180286 converter: one .h5ad file per sample
==============================================
Processes files in chunks with proper AnnData output.
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, vstack
import scanpy as sc
import gc

print("=" * 80)
print("GSE180286 CONVERTER - CHUNKED WITH PROPER ANNDATA OUTPUT")
print("=" * 80)

extract_dir = "GSE180286_extracted"
output_dir = "individual_samples"

os.makedirs(output_dir, exist_ok=True)

sample_files = sorted([f for f in os.listdir(extract_dir) 
                      if f.endswith('.txt.gz')])

print(f"\nFound {len(sample_files)} samples")
print("=" * 80)

processed_samples = []
failed_samples = []

def process_sample_chunked(filepath, output_file, sample_id, gsm_id):
    """Process a sample file in chunks"""
    
    print(f"    Reading in chunks of 500 cells...")
    
    # Read in chunks
    chunks = []
    cell_names = []
    gene_names = None
    
    chunk_iter = pd.read_csv(
        filepath,
        sep='\t',
        index_col=0,
        compression='gzip',
        chunksize=500,
        engine='c'
    )
    
    chunk_num = 0
    for chunk in chunk_iter:
        chunk_num += 1
        
        # Get gene names from first chunk
        if gene_names is None:
            gene_names = chunk.columns.tolist()
            print(f"    Found {len(gene_names):,} genes")
        
        print(f"      Chunk {chunk_num}: {len(chunk):,} cells...", end=' ')
        
        # Convert to sparse
        chunk_sparse = csr_matrix(chunk.values, dtype=np.float32)
        chunks.append(chunk_sparse)
        cell_names.extend(chunk.index.tolist())
        
        del chunk
        gc.collect()
        print("done")
    
    print(f"    Total: {len(cell_names):,} cells")
    print(f"    Combining {len(chunks)} chunks...")
    
    # Combine all chunks
    X_combined = vstack(chunks, format='csr')
    del chunks
    gc.collect()
    
    print(f"    Creating AnnData object...")
    
    # Create proper AnnData object
    adata = sc.AnnData(X=X_combined, dtype=np.float32)
    adata.var_names = gene_names
    adata.obs_names = cell_names
    adata.obs['sample_id'] = sample_id
    adata.obs['gsm_id'] = gsm_id
    
    # Add QC metrics
    print(f"    Calculating QC metrics...")
    adata.var['mt'] = [str(x).startswith('MT-') for x in adata.var_names]
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # Save using scanpy's writer
    print(f"    Saving with scanpy...")
    adata.write_h5ad(output_file, compression='gzip')
    
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    
    del adata
    del X_combined
    del cell_names
    del gene_names
    gc.collect()
    
    return n_cells, n_genes

# Process each sample
for i, file in enumerate(sample_files, 1):
    filepath = os.path.join(extract_dir, file)
    
    base_name = file.replace('.expression_matrix.txt.gz', '')
    parts = base_name.split('_')
    gsm_id = parts[0]
    sample_id = parts[1] if len(parts) > 1 else f'sample_{i}'
    
    output_file = os.path.join(output_dir, f"{gsm_id}_{sample_id}.h5ad")
    
    if os.path.exists(output_file):
        print(f"\n[{i}/{len(sample_files)}] {sample_id} - Already exists, skipping")
        processed_samples.append(output_file)
        continue
    
    print(f"\n[{i}/{len(sample_files)}] Processing: {sample_id}")
    print(f"    File: {file}")
    
    try:
        n_cells, n_genes = process_sample_chunked(filepath, output_file, sample_id, gsm_id)
        
        file_size = os.path.getsize(output_file) / (1024**2)
        print(f"    Complete: {file_size:.1f} MB")
        
        processed_samples.append(output_file)
        
    except MemoryError:
        print(f"    Out of memory.")
        failed_samples.append(file)
    except Exception as e:
        print(f"    Error: {e}")
        import traceback
        traceback.print_exc()
        failed_samples.append(file)

print(f"\n{'=' * 80}")
print(f"PROCESSING COMPLETE")
print(f"{'=' * 80}")
print(f"\nSuccessfully processed: {len(processed_samples)}/{len(sample_files)} samples")

if failed_samples:
    print(f"\nFailed samples ({len(failed_samples)}):")
    for f in failed_samples:
        print(f"  - {f}")

if len(processed_samples) == 0:
    print("\nNo samples were processed.")
    sys.exit(1)

print(f"\nFiles saved in: {output_dir}/")
print(f"\nNext step: python merge_samples.py")
print(f"\n{'=' * 80}\n")
