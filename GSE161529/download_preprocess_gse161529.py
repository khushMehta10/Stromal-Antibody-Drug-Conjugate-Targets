"""
DOWNLOAD AND PREPROCESS GSE161529
==================================
Complete pipeline to download and process breast cancer scRNA-seq data
Preserves all genes including stromal markers

Runtime: ~1-2 hours total
Memory: ~4-6 GB peak
"""

import subprocess
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import gzip
import gc
from scipy import sparse, io

# ============================================================================
# PART 1: DOWNLOAD DATA
# ============================================================================

def download_gse161529():
    """Download GSE161529 data from GEO using prefetch/fasterq-dump"""
    
    print("\n" + "="*80)
    print("DOWNLOADING GSE161529 FROM GEO")
    print("="*80)
    
    # Check if SRA toolkit is installed
    try:
        subprocess.run(['prefetch', '--version'], check=True, 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("SRA Toolkit found")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("\nSRA Toolkit not found.")
        print("\nInstall instructions:")
        print("  conda install -c bioconda sra-tools")
        print("  or visit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit")
        sys.exit(1)
    
    # Create download directory
    download_dir = Path('./geo_raw_data')
    download_dir.mkdir(exist_ok=True)
    
    print(f"\nDownload directory: {download_dir.absolute()}")
    
    # Get list of all SRR accessions for GSE161529
    # This dataset has 45 samples
    print("\nNOTE: This dataset is LARGE (~150GB compressed)")
    print("Download will take 1-2 hours depending on connection speed")
    
    response = input("\nProceed with download? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("Download canceled")
        sys.exit(0)
    
    # GSE161529 SRR accessions (first 10 samples as example)
    # Full list: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529
    srr_accessions = [
        'SRR13045891', 'SRR13045892', 'SRR13045893', 'SRR13045894', 'SRR13045895',
        'SRR13045896', 'SRR13045897', 'SRR13045898', 'SRR13045899', 'SRR13045900',
        # Add all 45 accessions here - truncated for brevity
        # See: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA682935
    ]
    
    print(f"\nWill download {len(srr_accessions)} samples")
    print("(Note: Full dataset has 45 samples - add remaining accessions to script)")
    
    for i, srr in enumerate(srr_accessions, 1):
        print(f"\n[{i}/{len(srr_accessions)}] Downloading {srr}...")
        
        # Prefetch
        cmd = ['prefetch', srr, '-O', str(download_dir)]
        try:
            subprocess.run(cmd, check=True)
            print(f"  Prefetched")
        except subprocess.CalledProcessError as e:
            print(f"  Failed: {e}")
            continue
        
        # Convert to FASTQ
        cmd = ['fasterq-dump', srr, '-O', str(download_dir), '--split-files']
        try:
            subprocess.run(cmd, check=True)
            print(f"  Converted to FASTQ")
        except subprocess.CalledProcessError as e:
            print(f"  Failed: {e}")
    
    print("\nDownload complete")
    return download_dir


def download_geo_series_matrix():
    """Alternative: Download pre-processed matrix files from GEO"""
    
    print("\n" + "="*80)
    print("DOWNLOADING GEO SERIES MATRIX (RECOMMENDED)")
    print("="*80)
    
    print("\nThis is faster and includes metadata.")
    print("You'll need to manually download from:")
    print("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529")
    print("\nLook for 'Supplementary file' links - download all MTX files")
    
    download_dir = Path('./geo_raw_data')
    download_dir.mkdir(exist_ok=True)
    
    print(f"\nDownload to: {download_dir.absolute()}")
    print("\nFiles needed:")
    print("  - GSE161529_*_matrix.mtx.gz (45 files)")
    print("  - GSE161529_*_barcodes.tsv.gz (45 files)")
    print("  - GSE161529_features.tsv.gz (1 shared file)")
    
    response = input("\nHave you downloaded these files? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("\nPlease download files first, then run this script again")
        sys.exit(0)
    
    return download_dir


# ============================================================================
# PART 2: PREPROCESS DATA (MINIMAL FILTERING)
# ============================================================================

def find_sample_files(data_dir):
    """Find all MTX sample files"""
    
    print("\n" + "="*80)
    print("DISCOVERING SAMPLES")
    print("="*80)
    
    data_path = Path(data_dir)
    
    # Find shared features file
    features_files = list(data_path.glob("*features.tsv.gz"))
    if not features_files:
        raise FileNotFoundError("No features.tsv.gz file found")
    
    shared_features = features_files[0]
    print(f"\nShared features: {shared_features.name}")
    
    # Find all matrix files
    matrix_files = sorted(data_path.glob("*matrix.mtx.gz"))
    print(f"Found {len(matrix_files)} matrix files")
    
    # Match matrices with barcodes
    samples = {}
    
    for matrix_file in matrix_files:
        # Get barcode file
        barcode_file = Path(str(matrix_file).replace('matrix.mtx.gz', 'barcodes.tsv.gz'))
        
        if barcode_file.exists():
            prefix = matrix_file.stem.replace('_matrix.mtx', '').replace('.mtx', '')
            samples[prefix] = {
                'matrix': matrix_file,
                'barcodes': barcode_file,
                'features': shared_features
            }
            print(f"  {prefix}: OK")
        else:
            print(f"  {matrix_file.name} - missing barcodes")
    
    if len(samples) == 0:
        raise FileNotFoundError("No complete sample sets found!")
    
    print(f"\nTotal samples: {len(samples)}")
    return samples


def load_mtx_sample(sample_prefix, files):
    """Load one MTX sample"""
    
    # Load matrix
    with gzip.open(files['matrix'], 'rb') as f:
        matrix = io.mmread(f)
    
    # Load barcodes and features
    barcodes = pd.read_csv(files['barcodes'], sep='\t', header=None)[0].values
    features_df = pd.read_csv(files['features'], sep='\t', header=None)
    gene_names = features_df[1].values if features_df.shape[1] >= 2 else features_df[0].values
    
    # Transpose if needed (MTX is usually genes x cells)
    if matrix.shape[0] == len(gene_names):
        X = sparse.csr_matrix(matrix.T, dtype=np.float32)
    else:
        X = sparse.csr_matrix(matrix, dtype=np.float32)
    
    # Create AnnData
    adata = sc.AnnData(X=X)
    adata.obs_names = [f"{sample_prefix}_{bc}" for bc in barcodes[:X.shape[0]]]
    adata.var_names = gene_names[:X.shape[1]]
    adata.var_names_make_unique()
    
    # Parse metadata from sample name
    # Format: GSM4909253_N-PM0092-Total
    parts = sample_prefix.split('_')
    gsm_id = parts[0] if len(parts) > 0 else sample_prefix
    sample_info = '_'.join(parts[1:]) if len(parts) > 1 else sample_prefix
    
    # Parse patient/sample info
    info_parts = sample_info.split('-')
    
    if len(info_parts) >= 2:
        subtype_code = info_parts[0]
        patient_id = info_parts[1]
        
        # Determine sample type
        if len(info_parts) >= 3:
            sample_suffix = info_parts[2].lower()
            sample_type = 'metastatic' if 'ln' in sample_suffix or 'met' in sample_suffix else 'primary'
        else:
            sample_type = 'primary'
        
        # Map cancer subtype
        if subtype_code == 'N':
            cancer_subtype = 'Normal'
        elif subtype_code in ['TN', 'B1']:
            cancer_subtype = 'TNBC'
        elif 'ER' in subtype_code:
            cancer_subtype = 'ER+'
        elif subtype_code == 'HER2':
            cancer_subtype = 'HER2+'
        else:
            cancer_subtype = 'Unknown'
    else:
        patient_id = sample_info
        sample_type = 'unknown'
        cancer_subtype = 'Unknown'
    
    # Add metadata
    adata.obs['sample_id'] = sample_prefix
    adata.obs['gsm_id'] = gsm_id
    adata.obs['patient_id'] = patient_id
    adata.obs['sample_type'] = sample_type
    adata.obs['cancer_subtype'] = cancer_subtype
    adata.obs['batch'] = f"batch_{cancer_subtype}"
    
    return adata


def load_all_samples(samples, batch_size=5):
    """Load all samples in batches"""
    
    print("\n" + "="*80)
    print("LOADING ALL SAMPLES")
    print("="*80)
    
    sample_list = list(samples.items())
    all_adatas = []
    
    for i in range(0, len(sample_list), batch_size):
        batch = sample_list[i:i+batch_size]
        print(f"\nBatch {i//batch_size + 1}/{(len(sample_list)-1)//batch_size + 1}")
        
        batch_adatas = []
        for prefix, files in batch:
            print(f"  Loading {prefix}...", end=' ')
            adata = load_mtx_sample(prefix, files)
            print(f"{adata.n_obs:,} cells")
            batch_adatas.append(adata)
            gc.collect()
        
        # Concatenate batch
        batch_combined = sc.concat(batch_adatas, join='outer', index_unique=None)
        all_adatas.append(batch_combined)
        
        del batch_adatas
        gc.collect()
    
    # Final concatenation
    print("\nConcatenating all batches...")
    adata_combined = sc.concat(all_adatas, join='outer', index_unique=None)
    
    del all_adatas
    gc.collect()
    
    print(f"\nCombined dataset:")
    print(f"  Cells: {adata_combined.n_obs:,}")
    print(f"  Genes: {adata_combined.n_vars:,}")
    print(f"  Patients: {adata_combined.obs['patient_id'].nunique()}")
    
    return adata_combined


def minimal_qc_filtering(adata):
    """MINIMAL filtering to preserve stromal markers"""
    
    print("\n" + "="*80)
    print("QUALITY CONTROL (MINIMAL FILTERING)")
    print("="*80)
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                               log1p=False, inplace=True)
    
    print("\nBefore filtering:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Median counts/cell: {adata.obs['total_counts'].median():.0f}")
    
    # MINIMAL cell filtering
    n_cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=500)
    adata = adata[adata.obs['pct_counts_mt'] <= 20, :].copy()
    
    print(f"\nRemoved {n_cells_before - adata.n_obs:,} low-quality cells")
    
    # Minimal gene filtering: keep genes detected in at least 10 cells
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=10)  # permissive threshold
    
    print(f"Removed {n_genes_before - adata.n_vars:,} genes (expressed in <10 cells)")
    
    # Verify stromal markers are present
    print("\n" + "="*80)
    print("VERIFYING STROMAL MARKERS")
    print("="*80)
    
    stromal_markers = ['FAP', 'PDPN', 'COL1A1', 'COL1A2', 'PDGFRA', 'PDGFRB', 
                       'DCN', 'LUM', 'VIM', 'ACTA2', 'POSTN', 'THY1']
    
    found = [g for g in stromal_markers if g in adata.var_names]
    missing = [g for g in stromal_markers if g not in adata.var_names]
    
    print(f"\nFound {len(found)}/{len(stromal_markers)} stromal markers:")
    for gene in found:
        print(f"  {gene}")
    
    if missing:
        print(f"\nMissing markers:")
        for gene in missing:
            print(f"  {gene}")
    
    print(f"\nAfter filtering:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    
    gc.collect()
    return adata


def normalize_and_cluster(adata):
    """Normalize and cluster"""
    
    print("\n" + "="*80)
    print("NORMALIZATION AND CLUSTERING")
    print("="*80)
    
    # Save raw counts
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize
    print("\nNormalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()
    
    # Find highly variable genes for dimensionality reduction
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch',
                                flavor='seurat_v3', layer='counts', subset=False)
    
    print(f"{adata.var['highly_variable'].sum()} highly variable genes")
    
    # PCA on HVGs
    print("\nRunning PCA...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.regress_out(adata_hvg, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
    
    # Transfer PCA back to full dataset
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
    adata.varm['PCs'] = np.zeros((adata.n_vars, 50))
    adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']
    adata.uns['pca'] = adata_hvg.uns['pca'].copy()
    
    del adata_hvg
    gc.collect()
    
    # Batch correction with Harmony (optional)
    try:
        print("\nBatch correction with Harmony...")
        import harmonypy as hm
        ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch', max_iter_harmony=20)
        adata.obsm['X_pca_harmony'] = ho.Z_corr.T
        use_rep = 'X_pca_harmony'
        print("Harmony complete")
    except ImportError:
        print("harmonypy not installed, skipping batch correction")
        adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()
        use_rep = 'X_pca'
    
    # Clustering
    print("\nClustering...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0, key_added='leiden')
    
    print(f"{adata.obs['leiden'].nunique()} Leiden clusters")
    
    # Basic cell type annotation (OPTIMIZED FOR LARGE DATASETS)
    print("\nAnnotating cell types (this may take 5-10 minutes)...")
    
    # Simple marker-based annotation
    markers = {
        'Epithelial': ['EPCAM', 'KRT19', 'KRT18', 'KRT8'],
        'Fibroblast': ['FAP', 'COL1A1', 'DCN', 'LUM'],
        'T_cell': ['CD3D', 'CD3E', 'CD3G'],
        'B_cell': ['CD79A', 'MS4A1', 'CD19'],
        'Myeloid': ['CD14', 'CD68', 'LYZ'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5']
    }
    
    cluster_annotations = {}
    clusters = sorted(adata.obs['leiden'].unique(), key=lambda x: int(x))
    
    for i, cluster in enumerate(clusters, 1):
        print(f"  [{i}/{len(clusters)}] Cluster {cluster}...", end=' ')
        
        cluster_mask = adata.obs['leiden'] == cluster
        n_cells = cluster_mask.sum()
        
        scores = {}
        
        for cell_type, genes in markers.items():
            available = [g for g in genes if g in adata.var_names]
            if available:
                # Get indices for faster subsetting
                cluster_indices = np.where(cluster_mask)[0]
                gene_indices = [adata.var_names.get_loc(g) for g in available]
                
                # Extract expression efficiently
                if hasattr(adata.X, 'toarray'):
                    # Sparse matrix
                    expr_sum = 0
                    for gene_idx in gene_indices:
                        expr_sum += adata.X[cluster_indices, gene_idx].toarray().mean()
                    scores[cell_type] = expr_sum / len(available)
                else:
                    # Dense matrix
                    expr = adata.X[np.ix_(cluster_indices, gene_indices)]
                    scores[cell_type] = expr.mean()
            else:
                scores[cell_type] = 0
        
        best_type = max(scores, key=scores.get)
        cluster_annotations[cluster] = best_type
        print(f"{best_type} ({n_cells:,} cells)")
        
        # Clear memory
        if i % 10 == 0:
            gc.collect()
    
    adata.obs['celltype_annotation'] = adata.obs['leiden'].map(cluster_annotations)
    
    print("\nCell type distribution:")
    for celltype, count in adata.obs['celltype_annotation'].value_counts().items():
        pct = 100 * count / adata.n_obs
        print(f"  {celltype:15s}: {count:7,} ({pct:4.1f}%)")
    
    gc.collect()


def save_dataset(adata, output_dir='./data'):
    """Save processed dataset"""
    
    print("\n" + "="*80)
    print("SAVING DATASET")
    print("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    output_file = output_path / 'dataset.h5ad'
    
    print(f"\nSaving to: {output_file}")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Size: ~{adata.n_obs * adata.n_vars * 4 / 1e9:.2f} GB")
    
    adata.write_h5ad(output_file, compression='gzip')
    
    print(f"\nDataset saved.")
    
    # Save summary
    summary = {
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'n_patients': int(adata.obs['patient_id'].nunique()),
        'n_clusters': int(adata.obs['leiden'].nunique()),
        'sample_types': adata.obs['sample_type'].value_counts().to_dict(),
        'cell_types': adata.obs['celltype_annotation'].value_counts().to_dict()
    }
    
    import json
    with open(output_path / 'dataset_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved: {output_path / 'dataset_summary.json'}")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Main pipeline"""
    
    print("\n" + "="*80)
    print("GSE161529 DOWNLOAD AND PREPROCESSING PIPELINE")
    print("="*80)
    print("\nThis script will:")
    print("  1. Download GSE161529 data from GEO")
    print("  2. Process with MINIMAL gene filtering")
    print("  3. Preserve all stromal marker genes")
    print("  4. Create dataset.h5ad for analysis")
    
    # Step 1: Download
    print("\n" + "="*80)
    print("STEP 1: DOWNLOAD DATA")
    print("="*80)
    print("\nChoose download method:")
    print("  1. Manual download (RECOMMENDED)")
    print("  2. Automatic download (requires SRA toolkit)")
    
    choice = input("\nEnter choice (1 or 2): ").strip()
    
    if choice == '1':
        data_dir = download_geo_series_matrix()
    elif choice == '2':
        data_dir = download_gse161529()
    else:
        print("Invalid choice")
        sys.exit(1)
    
    # Step 2: Load samples
    samples = find_sample_files(data_dir)
    adata = load_all_samples(samples, batch_size=5)
    
    # Step 3: QC and filtering
    adata = minimal_qc_filtering(adata)
    
    # Step 4: Normalize and cluster
    normalize_and_cluster(adata)
    
    # Step 5: Save
    save_dataset(adata)
    
    print("\n" + "="*80)
    print("PIPELINE COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("  1. Verify stromal markers: python check_gse161529.py")
    print("  2. Run analysis pipeline: python run_complete_pipeline.py")
    print("\n" + "="*80)


if __name__ == "__main__":
    main()
