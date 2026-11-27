#!/usr/bin/env python3
"""
Merge Individual Sample Files
==============================
Merges all individual .h5ad files into a single dataset.
Uses memory-efficient batched merging for large datasets.
"""

import os
import sys
import scanpy as sc
import gc

print("=" * 80)
print("MERGING INDIVIDUAL SAMPLES")
print("=" * 80)

output_dir = "individual_samples"

# Check if directory exists
if not os.path.exists(output_dir):
    print(f"\nError: Directory '{output_dir}' not found.")
    print("Run convert_samples.py first.")
    sys.exit(1)

# Get all sample files
sample_files = sorted([f for f in os.listdir(output_dir) if f.endswith('.h5ad')])

if len(sample_files) == 0:
    print(f"\nError: No .h5ad files found in '{output_dir}'.")
    sys.exit(1)

print(f"\nFound {len(sample_files)} sample files to merge")
print(f"Directory: {output_dir}/")

# Show sample files
print(f"\nSample files:")
for i, f in enumerate(sample_files[:10], 1):
    size = os.path.getsize(os.path.join(output_dir, f)) / (1024**2)
    print(f"  {i:2d}. {f:50s} ({size:6.1f} MB)")
if len(sample_files) > 10:
    print(f"  ... and {len(sample_files) - 10} more")

print("\n" + "=" * 80)

# Estimate memory requirements
try:
    import psutil
    available_mem = psutil.virtual_memory().available / (1024**3)
    print(f"\nAvailable RAM: {available_mem:.1f} GB")
except:
    print("\nNote: Install psutil to see memory info (pip install psutil)")
    available_mem = None

# Decide merging strategy
if len(sample_files) <= 5:
    # Small number of samples - merge all at once
    print(f"\nStrategy: Merge all {len(sample_files)} samples at once")
    print("=" * 80)
    
    samples = []
    for i, f in enumerate(sample_files, 1):
        filepath = os.path.join(output_dir, f)
        print(f"\n[{i}/{len(sample_files)}] Loading {f}...")
        adata = sc.read_h5ad(filepath)
        print(f"    Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        samples.append(adata)
    
    print(f"\nMerging {len(samples)} samples...")
    adata_combined = sc.concat(samples, join='outer', index_unique='-')
    
else:
    # Many samples - merge in batches
    batch_size = 3
    print(f"\nStrategy: Merge in batches of {batch_size} (memory-efficient)")
    print("=" * 80)
    
    # First pass: merge into batches
    batch_files = []
    n_batches = (len(sample_files) + batch_size - 1) // batch_size
    
    print(f"\nStep 1: Creating {n_batches} batch files...")
    
    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, len(sample_files))
        batch = sample_files[start_idx:end_idx]
        
        print(f"\n  Batch {batch_idx + 1}/{n_batches}: Merging {len(batch)} samples...")
        
        batch_samples = []
        for f in batch:
            filepath = os.path.join(output_dir, f)
            print(f"    Loading {f}...")
            adata = sc.read_h5ad(filepath)
            batch_samples.append(adata)
        
        print(f"    Concatenating...")
        batch_merged = sc.concat(batch_samples, join='outer', index_unique='-')
        
        # Save batch
        batch_file = f"temp_batch_{batch_idx + 1}.h5ad"
        print(f"    Saving batch to {batch_file}...")
        batch_merged.write(batch_file, compression='gzip')
        batch_files.append(batch_file)
        
        print(f"    Batch {batch_idx + 1} complete: {batch_merged.n_obs:,} cells")
        
        # Clean up
        del batch_samples
        del batch_merged
        gc.collect()
    
    # Second pass: merge all batches
    print(f"\nStep 2: Merging {len(batch_files)} batches into final dataset...")
    
    final_samples = []
    for i, batch_file in enumerate(batch_files, 1):
        print(f"  [{i}/{len(batch_files)}] Loading {batch_file}...")
        adata = sc.read_h5ad(batch_file)
        print(f"      {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        final_samples.append(adata)
    
    print(f"\n  Merging batches...")
    adata_combined = sc.concat(final_samples, join='outer', index_unique='-')
    
    # Clean up batch files
    print(f"\n  Cleaning up temporary batch files...")
    for batch_file in batch_files:
        os.remove(batch_file)
        print(f"    Removed {batch_file}")

# Save final merged file
print("\n" + "=" * 80)
print("SAVING FINAL DATASET")
print("=" * 80)

final_output = "GSE180286_breast_cancer_raw.h5ad"

print(f"\nFinal dataset:")
print(f"  Cells:   {adata_combined.n_obs:,}")
print(f"  Genes:   {adata_combined.n_vars:,}")
print(f"  Samples: {adata_combined.obs['sample_id'].nunique()}")

print(f"\nSaving to: {final_output}")
print("(This may take 5-10 minutes for large datasets...)")

adata_combined.write(final_output, compression='gzip')

file_size = os.path.getsize(final_output) / (1024**3)
print(f"Saved: {file_size:.2f} GB")

# Summary
print("\n" + "=" * 80)
print("MERGE COMPLETE")
print("=" * 80)

print(f"\nFinal file: {final_output}")
print(f"File size:  {file_size:.2f} GB")
print(f"\nDataset summary:")
print(f"  • Total cells:  {adata_combined.n_obs:,}")
print(f"  • Total genes:  {adata_combined.n_vars:,}")
print(f"  • Total samples: {adata_combined.obs['sample_id'].nunique()}")

print(f"\nSample breakdown:")
sample_counts = adata_combined.obs['sample_id'].value_counts()
for sample in sample_counts.index[:10]:
    print(f"  {sample:30s}: {sample_counts[sample]:8,} cells")
if len(sample_counts) > 10:
    print(f"  ... and {len(sample_counts) - 10} more samples")

print(f"\nReady for analysis.")
print(f"\nNext steps:")
print(f"  import scanpy as sc")
print(f"  adata = sc.read_h5ad('{final_output}')")

print("\n" + "=" * 80 + "\n")
