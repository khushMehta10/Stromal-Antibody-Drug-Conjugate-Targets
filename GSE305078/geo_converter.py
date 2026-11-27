"""
Convert GSE305078 data to h5ad format (sparse merge)
====================================================
Uses optimized sparse concatenation without reindexing each matrix
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import scipy.io
import scipy.sparse as sp
import re
from collections import defaultdict
import gc
import anndata as ad

def find_sample_groups(data_dir):
    """Group files by sample ID"""
    data_path = Path(data_dir)
    all_files = list(data_path.glob('*.gz'))
    
    if not all_files:
        return {}
    
    sample_groups = defaultdict(dict)
    
    for filepath in all_files:
        filename = filepath.name
        match = re.match(r'(GSM\d+_[EA]\d+)_(barcodes|features|matrix)', filename)
        if match:
            sample_id = match.group(1)
            file_type = match.group(2)
            sample_groups[sample_id][file_type] = filepath
    
    return sample_groups


def load_10x_sparse(sample_id, file_dict):
    """Load sample and return sparse matrix + metadata"""
    required = ['matrix', 'barcodes', 'features']
    if not all(f in file_dict for f in required):
        return None
    
    try:
        matrix_file = file_dict['matrix']
        barcode_file = file_dict['barcodes']
        feature_file = file_dict['features']
        
        # Load matrix
        with gzip.open(matrix_file, 'rb') as f:
            mat = scipy.io.mmread(f).T.tocsr()
        
        # Load barcodes
        barcodes = pd.read_csv(barcode_file, header=None, compression='gzip')[0].values
        
        # Load features
        features = pd.read_csv(feature_file, header=None, sep='\t', compression='gzip')
        gene_names = features[1].values if features.shape[1] >= 2 else features[0].values
        
        # Add sample prefix to barcodes
        barcodes = [f"{sample_id}_{bc}" for bc in barcodes]
        
        # Extract patient ID
        patient_match = re.search(r'([EA]\d+)', sample_id)
        patient_id = patient_match.group(1) if patient_match else sample_id
        
        return {
            'matrix': mat,
            'barcodes': barcodes,
            'genes': gene_names,
            'sample_id': sample_id,
            'patient_id': patient_id,
            'n_cells': mat.shape[0]
        }
        
    except Exception as e:
        print(f"    Error: {e}")
        return None


def merge_samples_in_batches(sample_data_list, batch_size=4):
    """Merge samples in small batches using scanpy's concat"""
    print("\nMerging samples in batches...")
    print(f"   Batch size: {batch_size} samples\n")
    
    batch_adatas = []
    
    # Process in batches
    for i in range(0, len(sample_data_list), batch_size):
        batch = sample_data_list[i:i+batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(sample_data_list) + batch_size - 1) // batch_size
        
        print(f"  Batch {batch_num}/{total_batches}: Processing {len(batch)} samples...")
        
        # Create AnnData objects for this batch
        adatas = []
        for data in batch:
            obs = pd.DataFrame({
                'sample_id': data['sample_id'],
                'patient_id': data['patient_id']
            }, index=data['barcodes'])
            
            var = pd.DataFrame(index=data['genes'])
            
            adata = ad.AnnData(X=data['matrix'], obs=obs, var=var)
            adata.var_names_make_unique()
            adatas.append(adata)
            
            # Clean up
            del data['matrix']
        
        # Concatenate batch
        if len(adatas) == 1:
            batch_adata = adatas[0]
        else:
            batch_adata = sc.concat(adatas, axis=0, join='outer', merge='unique')
        
        print(f"    {batch_adata.n_obs:,} cells × {batch_adata.n_vars:,} genes")
        
        batch_adatas.append(batch_adata)
        
        # Clean up
        del adatas
        gc.collect()
    
    # Now merge all batches
    print(f"\n  Final merge of {len(batch_adatas)} batches...")
    
    if len(batch_adatas) == 1:
        return batch_adatas[0]
    
    # Merge batches pairwise to control memory
    while len(batch_adatas) > 1:
        print(f"    Merging {len(batch_adatas)} intermediate files...")
        next_level = []
        
        for i in range(0, len(batch_adatas), 2):
            if i + 1 < len(batch_adatas):
                # Merge pair
                merged = sc.concat([batch_adatas[i], batch_adatas[i+1]], 
                                 axis=0, join='outer', merge='unique')
                next_level.append(merged)
                gc.collect()
            else:
                # Odd one out
                next_level.append(batch_adatas[i])
        
        # Clear old batch_adatas
        del batch_adatas
        batch_adatas = next_level
        gc.collect()
        
        print(f"    → {len(batch_adatas)} files remaining")
    
    return batch_adatas[0]


def load_and_merge_all(data_dir):
    """Load all samples and merge efficiently"""
    print("="*80)
    print("LOADING ALL SAMPLES")
    print("="*80)
    
    sample_groups = find_sample_groups(data_dir)
    if not sample_groups:
        return None
    
    sample_ids = sorted(sample_groups.keys())
    print(f"\nFound {len(sample_ids)} samples\n")
    
    # Load all samples as sparse data
    sample_data_list = []
    
    for i, sample_id in enumerate(sample_ids, 1):
        print(f"[{i:2d}/{len(sample_ids)}] Loading {sample_id}...", end=" ")
        
        file_dict = sample_groups[sample_id]
        data = load_10x_sparse(sample_id, file_dict)
        
        if data is None:
            print("FAILED")
            continue
        
        print(f"{data['n_cells']:,} cells")
        sample_data_list.append(data)
    
    if not sample_data_list:
        return None
    
    print(f"\nLoaded {len(sample_data_list)} samples")
    
    # Merge in batches (smaller batch size for lower memory usage)
    adata = merge_samples_in_batches(sample_data_list, batch_size=2)
    
    # Clean up
    del sample_data_list
    gc.collect()
    
    print(f"\nFinal dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    return adata


def basic_qc_and_cluster(adata):
    """QC and clustering"""
    print("\n" + "="*80)
    print("QC AND CLUSTERING")
    print("="*80)
    
    # QC
    print("\nQC metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    if adata.var['mt'].sum() > 0:
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        print(f"  {adata.var['mt'].sum()} mitochondrial genes")
    
    print(f"\nFiltering...")
    print(f"  Before: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Sample types
    if 'patient_id' in adata.obs.columns:
        adata.obs['sample_type'] = adata.obs['patient_id'].apply(
            lambda x: 'early_stage' if str(x).startswith('E') else 
                     'advanced' if str(x).startswith('A') else 'unknown'
        )
        
        print(f"\nSample types:")
        for stype, count in adata.obs['sample_type'].value_counts().items():
            pct = 100 * count / len(adata)
            print(f"  {stype}: {count:,} ({pct:.1f}%)")
    
    gc.collect()
    
    # Clustering
    print("\nClustering...")
    print("  Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print("  Variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
    
    # Subset to HVGs to reduce memory for PCA
    print("  Subsetting to HVGs for dimensional reduction...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    
    print("  Scaling HVGs only...")
    sc.pp.scale(adata_hvg, max_value=10)
    
    print("  PCA...")
    sc.tl.pca(adata_hvg, n_comps=50)
    
    # Transfer PCA back to full dataset
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca'] = adata_hvg.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.n_vars, 50))
    adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']
    
    del adata_hvg
    gc.collect()
    
    print("  Neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    
    print("  UMAP...")
    sc.tl.umap(adata)
    
    print("  Leiden...")
    sc.tl.leiden(adata, resolution=0.5)
    
    print(f"{adata.obs['leiden'].nunique()} clusters")
    
    gc.collect()
    return adata


def main():
    """Main workflow"""
    print("\n" + "="*80)
    print("GSE305078 CONVERSION (BATCH MERGE)")
    print("="*80)
    print("\nMerging in small batches with pairwise final merge")
    
    data_dir = './GSE305078_data'
    output_dir = Path('./data')
    output_dir.mkdir(exist_ok=True)
    
    if not Path(data_dir).exists():
        print(f"\nData directory not found: {data_dir}")
        return
    
    # Load and merge all samples
    adata = load_and_merge_all(data_dir)
    
    if adata is None:
        print("\nFailed to load data")
        return
    
    # Save raw combined
    raw_file = output_dir / 'combined_raw.h5ad'
    print(f"\nSaving raw combined data: {raw_file}")
    adata.write_h5ad(raw_file)
    print(f"Saved ({raw_file.stat().st_size / 1e9:.2f} GB)")
    
    # QC and clustering
    adata = basic_qc_and_cluster(adata)
    
    # Save final
    final_file = output_dir / 'GSE305078_dataset.h5ad'
    print(f"\nSaving final dataset: {final_file}")
    adata.write_h5ad(final_file)
    
    print(f"Saved ({final_file.stat().st_size / 1e9:.2f} GB)")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")
    print(f"Patients: {adata.obs['patient_id'].nunique()}")
    print(f"Clusters: {adata.obs['leiden'].nunique()}")
    
    if 'sample_type' in adata.obs.columns:
        print(f"\nSample types:")
        for stype, count in adata.obs['sample_type'].value_counts().items():
            pct = 100 * count / len(adata)
            print(f"  {stype}: {count:,} ({pct:.1f}%)")
    
    print("\n" + "="*80)
    print("Conversion complete.")
    print("="*80)
    print("\nNext: python run_GSE305078.py")


if __name__ == "__main__":
    main()
