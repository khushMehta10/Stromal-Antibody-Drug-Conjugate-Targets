"""
Memory-efficient preprocessing for GSE180286
============================================
Processes 15 raw expression matrices from GSE180286 in chunks
Uses direct sparse matrix construction to avoid dense intermediate arrays

Runtime: ~45-60 minutes
Memory: ~1-2 GB (reduced from 6+ GB)
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import gc
from scipy import sparse

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white')

def load_sample_sparse(filepath):
    """
    Load one expression matrix file as sparse matrix (memory efficient)
    Constructs sparse matrix directly without creating dense array
    """
    print(f"\n  Loading: {filepath.name}")
    
    # Read the gzipped file
    with gzip.open(filepath, 'rt') as f:
        # Read header first
        header = f.readline().strip().split('\t')
        cell_barcodes = header[1:]  # Skip first column (gene names column label)
        n_cells = len(cell_barcodes)
        
        # Lists for sparse matrix construction
        rows = []  # cell indices
        cols = []  # gene indices
        data = []  # expression values
        gene_names = []
        
        print(f"    Expected: {n_cells} cells")
        
        for gene_idx, line in enumerate(f):
            # Skip empty lines
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')
            gene_name = parts[0]
            gene_names.append(gene_name)
            
            # Handle column mismatch
            n_values = min(len(parts) - 1, n_cells)
            if len(parts) - 1 != n_cells:
                print(f"\n    Warning: Gene {gene_name} has {len(parts)-1} values, using first {n_values}")
            
            # Add non-zero values to sparse matrix
            for cell_idx in range(n_values):
                val_str = parts[cell_idx + 1]
                if val_str and val_str != '0':
                    try:
                        val = float(val_str)
                        if val != 0:
                            rows.append(cell_idx)
                            cols.append(gene_idx)
                            data.append(val)
                    except ValueError:
                        pass  # Skip invalid values
            
            # Progress indicator
            if (gene_idx + 1) % 5000 == 0:
                print(f"    Read {gene_idx+1} genes, {len(data):,} non-zero values...", end='\r')
        
        n_genes = len(gene_names)
        print(f"    Shape: {n_genes} genes x {n_cells} cells                    ")
        print(f"    Non-zero values: {len(data):,} ({100*len(data)/(n_genes*n_cells):.2f}% dense)")
    
    # Create sparse matrix directly (cells x genes format)
    print(f"    Building sparse matrix...")
    X_sparse = sparse.coo_matrix(
        (data, (rows, cols)), 
        shape=(n_cells, n_genes),
        dtype=np.float32
    ).tocsr()
    
    # Clear intermediate data
    del rows, cols, data
    gc.collect()
    
    print(f"    Sparse matrix: {X_sparse.data.nbytes / 1e6:.1f} MB")
    
    # Create AnnData object with sparse matrix
    adata = sc.AnnData(X=X_sparse)
    adata.var_names = gene_names
    adata.obs_names = cell_barcodes
    
    # Extract sample information from filename
    sample_id = filepath.stem.replace('.expression_matrix.txt', '')
    parts = sample_id.split('_')
    
    if len(parts) >= 2:
        gsm_id = parts[0]
        patient_sample = parts[1]
        
        if '-' in patient_sample:
            patient_id, replicate = patient_sample.rsplit('-', 1)
        else:
            patient_id = patient_sample
            replicate = '1'
        
        year = patient_id[-4:]
        batch = f"batch_{year}"
        
        if patient_id[0] in ['A', 'B', 'C']: #### WRONG
            sample_type = 'primary'
        else:
            sample_type = 'metastatic'
    else:
        gsm_id = sample_id
        patient_id = sample_id
        replicate = '1'
        batch = 'batch_unknown'
        sample_type = 'unknown'
    
    # Add metadata
    adata.obs['sample_id'] = sample_id
    adata.obs['gsm_id'] = gsm_id
    adata.obs['patient_id'] = patient_id
    adata.obs['replicate'] = replicate
    adata.obs['batch'] = batch
    adata.obs['sample_type'] = sample_type
    
    print(f"    Patient: {patient_id}, Type: {sample_type}, Batch: {batch}")
    
    return adata


def load_samples_in_batches(data_dir, batch_size=5):
    """
    Load samples in smaller batches, process, then merge
    """
    print("\n" + "="*80)
    print("LOADING SAMPLES IN BATCHES (memory-efficient)")
    print("="*80)
    
    # Find all expression matrix files
    pattern = "*expression_matrix.txt.gz"
    files = sorted(Path(data_dir).glob(pattern))
    
    print(f"\nFound {len(files)} sample files")
    print(f"Processing in batches of {batch_size} samples")
    
    if len(files) == 0:
        raise FileNotFoundError(f"No files matching {pattern} in {data_dir}")
    
    # Process in batches
    all_batches = []
    
    for i in range(0, len(files), batch_size):
        batch_files = files[i:i+batch_size]
        print(f"\n--- Batch {i//batch_size + 1}/{(len(files)-1)//batch_size + 1} ---")
        
        # Load batch
        adatas = []
        for filepath in batch_files:
            adata = load_sample_sparse(filepath)
            adatas.append(adata)
            gc.collect()
        
        # Concatenate batch
        print(f"\n  Concatenating {len(adatas)} samples...")
        batch_combined = sc.concat(adatas, join='outer', index_unique='_')
        
        # Clear individual samples
        del adatas
        gc.collect()
        
        print(f"  Batch shape: {batch_combined.n_obs} cells x {batch_combined.n_vars} genes")
        print(f"  Batch memory: {batch_combined.X.data.nbytes / 1e6:.1f} MB")
        
        all_batches.append(batch_combined)
    
    # Merge all batches
    print("\n" + "="*80)
    print("MERGING ALL BATCHES")
    print("="*80)
    
    adata_combined = sc.concat(all_batches, join='outer', index_unique='_')
    
    # Clear batches
    del all_batches
    gc.collect()
    
    print(f"\nCombined dataset:")
    print(f"  Total cells: {adata_combined.n_obs:,}")
    print(f"  Total genes: {adata_combined.n_vars:,}")
    print(f"  Patients: {adata_combined.obs['patient_id'].nunique()}")
    print(f"  Memory: {adata_combined.X.data.nbytes / 1e6:.1f} MB (sparse)")
    
    # Sample type distribution
    print(f"\nSample types:")
    for sample_type, count in adata_combined.obs['sample_type'].value_counts().items():
        pct = 100 * count / adata_combined.n_obs
        print(f"  {sample_type:15s}: {count:7,} cells ({pct:5.1f}%)")
    
    return adata_combined


def calculate_qc_metrics(adata):
    """Calculate QC metrics"""
    print("\n" + "="*80)
    print("CALCULATING QC METRICS")
    print("="*80)
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    print("\nQC distributions:")
    print(f"  Genes per cell: median={adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Counts per cell: median={adata.obs['total_counts'].median():.0f}")
    print(f"  MT percentage: median={adata.obs['pct_counts_mt'].median():.1f}%")


def filter_cells_and_genes(adata):
    """Filter low-quality cells and genes"""
    print("\n" + "="*80)
    print("FILTERING")
    print("="*80)
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells
    print("\nFiltering cells...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=500)
    adata = adata[adata.obs['total_counts'] <= 50000, :].copy()
    adata = adata[adata.obs['pct_counts_mt'] <= 20, :].copy()
    
    print(f"  Removed {n_cells_before - adata.n_obs:,} cells")
    
    # Filter genes
    print("\nFiltering genes...")
    sc.pp.filter_genes(adata, min_cells=3)
    
    print(f"  Removed {n_genes_before - adata.n_vars:,} genes")
    print(f"\nFinal: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    gc.collect()
    return adata


def normalize_and_process(adata):
    """Normalize and find variable genes"""
    print("\n" + "="*80)
    print("NORMALIZATION & VARIABLE GENES")
    print("="*80)
    
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize
    print("\nNormalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()
    
    # Find variable genes
    print("\nFinding variable genes...")
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=2000,
        batch_key='batch',
        flavor='seurat_v3',
        layer='counts',
        subset=False  # Don't subset yet, keep for later
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"{n_hvg} highly variable genes")
    
    gc.collect()


def scale_pca_batch_correct(adata):
    """Scale, PCA, and batch correction on HVG only"""
    print("\n" + "="*80)
    print("SCALING, PCA & BATCH CORRECTION")
    print("="*80)
    
    # Subset to HVG for processing (saves memory)
    print("\nSubsetting to highly variable genes...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    print(f"  Working with {adata_hvg.n_vars} genes")
    
    # Regress and scale
    print("\nRegressing out technical variation...")
    sc.pp.regress_out(adata_hvg, ['total_counts', 'pct_counts_mt'])
    
    print("Scaling...")
    sc.pp.scale(adata_hvg, max_value=10)
    
    # PCA
    print("Running PCA...")
    sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
    
    # Copy PCA back to full dataset
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
    adata.varm['PCs'] = np.zeros((adata.n_vars, 50))
    adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']
    adata.uns['pca'] = adata_hvg.uns['pca'].copy()
    
    # Batch correction with Harmony
    print("\nBatch correction...")
    try:
        import harmonypy as hm
        ho = hm.run_harmony(
            adata.obsm['X_pca'], 
            adata.obs, 
            'batch',
            max_iter_harmony=20
        )
        adata.obsm['X_pca_harmony'] = ho.Z_corr.T
        print("Harmony complete")
    except ImportError:
        print("harmonypy not installed - skipping")
        adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()
    
    del adata_hvg
    gc.collect()


def cluster_and_annotate(adata):
    """Clustering, UMAP, and basic annotation"""
    print("\n" + "="*80)
    print("CLUSTERING & ANNOTATION")
    print("="*80)
    
    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    
    print("\nComputing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep=use_rep)
    
    print("Computing UMAP...")
    sc.tl.umap(adata)
    
    print("Leiden clustering...")
    sc.tl.leiden(adata, resolution=1.0, key_added='leiden')
    
    n_clusters = adata.obs['leiden'].nunique()
    print(f"{n_clusters} clusters")
    
    # Calculate marker genes
    print("\nFinding marker genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Simple annotation
    markers = {
        'Epithelial': ['EPCAM', 'KRT19', 'KRT18'],
        'Fibroblast': ['FAP', 'COL1A1', 'DCN'],
        'T_cell': ['CD3D', 'CD3E'],
        'B_cell': ['CD79A', 'MS4A1'],
        'Myeloid': ['CD14', 'CD68'],
    }
    
    cluster_annotations = {}
    for cluster in adata.obs['leiden'].unique():
        cluster_cells = adata[adata.obs['leiden'] == cluster]
        scores = {}
        for cell_type, genes in markers.items():
            available_genes = [g for g in genes if g in adata.var_names]
            if available_genes:
                scores[cell_type] = cluster_cells[:, available_genes].X.mean()
            else:
                scores[cell_type] = 0
        cluster_annotations[cluster] = max(scores, key=scores.get)
    
    adata.obs['celltype_annotation'] = adata.obs['leiden'].map(cluster_annotations)
    adata.obs['major_celltype'] = adata.obs['leiden'].astype(int)
    
    print("\nCell type distribution:")
    for celltype, count in adata.obs['celltype_annotation'].value_counts().items():
        pct = 100 * count / adata.n_obs
        print(f"  {celltype:15s}: {count:7,} ({pct:4.1f}%)")


def save_dataset(adata, output_dir):
    """Save processed dataset"""
    print("\n" + "="*80)
    print("SAVING")
    print("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    output_file = output_path / 'dataset.h5ad'
    
    print(f"\nSaving to: {output_file}")
    adata.write_h5ad(output_file, compression='gzip')
    
    print(f"\nComplete.")
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")


def main():
    """Main pipeline"""
    print("\n" + "="*80)
    print("MEMORY-EFFICIENT PREPROCESSING PIPELINE")
    print("="*80)
    
    DATA_DIR = './geo_raw_data'
    OUTPUT_DIR = './data'
    
    # Load in batches (reduced from 5 to 3 for safety)
    adata = load_samples_in_batches(DATA_DIR, batch_size=2)
    
    # QC metrics
    calculate_qc_metrics(adata)
    
    # Filter
    adata = filter_cells_and_genes(adata)
    
    # Normalize
    normalize_and_process(adata)
    
    # Scale, PCA, batch correction
    scale_pca_batch_correct(adata)
    
    # Cluster and annotate
    cluster_and_annotate(adata)
    
    # Save
    save_dataset(adata, OUTPUT_DIR)
    
    print("\n" + "="*80)
    print("COMPLETE")
    print("="*80)
    print("\nNext: python run_complete_pipeline.py")


if __name__ == "__main__":
    main()
