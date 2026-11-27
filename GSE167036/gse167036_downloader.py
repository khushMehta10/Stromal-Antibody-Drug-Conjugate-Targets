"""
GSE167036 Download and Preparation for Stromal ADC Target Discovery
====================================================================
Downloads and prepares GSE167036: Paired primary and lymph node metastatic breast cancer
118,845 cells from 8 patients with paired primary and metastatic samples

Files needed:
- GSE167036_cell_counts_matrix.mtx.gz (expression matrix)
- GSE167036_features.csv.gz (gene names)
- GSE167036_barcodes.csv.gz (cell barcodes)
- GSE167036_meta_all.csv.gz (cell metadata)
"""

import subprocess
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import gzip
from scipy.io import mmread
from scipy.sparse import csr_matrix

def check_dependencies():
    """Check required packages"""
    print("="*80)
    print("CHECKING DEPENDENCIES")
    print("="*80)
    
    required = {
        'scanpy': 'scanpy',
        'pandas': 'pandas', 
        'numpy': 'numpy',
        'scipy': 'scipy',
        'matplotlib': 'matplotlib',
        'seaborn': 'seaborn'
    }
    
    missing = []
    for package, name in required.items():
        try:
            __import__(package)
            print(f"  {name}: OK")
        except ImportError:
            print(f"  {name} - MISSING")
            missing.append(name)
    
    if missing:
        print(f"\nInstall missing packages:")
        print(f"  pip install {' '.join(missing)}")
        return False
    
    print("\nAll dependencies available")
    return True


def download_geo_files(gse_id="GSE167036"):
    """
    Download specific files for GSE167036
    """
    print("\n" + "="*80)
    print(f"DOWNLOADING {gse_id} FILES")
    print("="*80)
    
    # Create download directory
    download_dir = Path(f'./{gse_id}_download')
    download_dir.mkdir(exist_ok=True)
    
    print(f"\nDownload directory: {download_dir.absolute()}")
    
    # Required files for GSE167036
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167036/suppl/"
    
    required_files = [
        'GSE167036_cell_counts_matrix.mtx.gz',  # 555 MB - expression data
        'GSE167036_features.csv.gz',             # 160 KB - gene names
        'GSE167036_barcodes.csv.gz',             # 753 KB - cell IDs
        'GSE167036_meta_all.csv.gz'              # 3.8 MB - cell metadata
    ]
    
    print("\nDownloading required files...")
    
    all_success = True
    for filename in required_files:
        output_file = download_dir / filename
        
        # Skip if already downloaded
        if output_file.exists() and output_file.stat().st_size > 1000:
            print(f"  {filename} (already exists)")
            continue
        
        file_url = base_url + filename
        print(f"  Downloading {filename}...", end=' ')
        
        try:
            cmd = ['wget', '-q', '-O', str(output_file), file_url]
            result = subprocess.run(cmd, timeout=600, capture_output=True)
            
            if result.returncode == 0 and output_file.exists() and output_file.stat().st_size > 1000:
                size_mb = output_file.stat().st_size / 1e6
                print(f"({size_mb:.1f} MB)")
            else:
                print(f"FAILED")
                all_success = False
        except FileNotFoundError:
            print("wget not found.")
            all_success = False
            break
        except subprocess.TimeoutExpired:
            print("timeout")
            all_success = False
        except Exception as e:
            print(f"failed: {e}")
            all_success = False
    
    if not all_success:
        print("\n" + "="*80)
        print("MANUAL DOWNLOAD REQUIRED")
        print("="*80)
        print(f"""
Download files manually from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167036

Required files:
  1. GSE167036_cell_counts_matrix.mtx.gz (555 MB)
  2. GSE167036_features.csv.gz (160 KB)
  3. GSE167036_barcodes.csv.gz (753 KB)
  4. GSE167036_meta_all.csv.gz (3.8 MB)

Save to: {download_dir.absolute()}

Then run this script again.
""")
        return None
    
    print("\nAll files downloaded.")
    return download_dir


def load_mtx_data(download_dir):
    """
    Load Matrix Market format data
    """
    print("\n" + "="*80)
    print("LOADING EXPRESSION DATA")
    print("="*80)
    
    # File paths
    mtx_file = download_dir / 'GSE167036_cell_counts_matrix.mtx.gz'
    features_file = download_dir / 'GSE167036_features.csv.gz'
    barcodes_file = download_dir / 'GSE167036_barcodes.csv.gz'
    
    # Check files exist
    for f in [mtx_file, features_file, barcodes_file]:
        if not f.exists():
            raise FileNotFoundError(f"Missing file: {f.name}")
    
    print("\n1. Loading gene names...")
    features_df = pd.read_csv(features_file)
    print(f"   {len(features_df):,} genes")
    print(f"   Columns: {features_df.columns.tolist()}")
    
    # Gene names might be in different columns
    if 'gene_name' in features_df.columns:
        gene_names = features_df['gene_name'].tolist()
    elif 'feature_name' in features_df.columns:
        gene_names = features_df['feature_name'].tolist()
    elif len(features_df.columns) >= 2:
        gene_names = features_df.iloc[:, 1].tolist()  # Usually 2nd column
    else:
        gene_names = features_df.iloc[:, 0].tolist()  # Fallback to first column
    
    print(f"   First 5 genes: {gene_names[:5]}")
    
    print("\n2. Loading cell barcodes...")
    barcodes_df = pd.read_csv(barcodes_file)
    print(f"   {len(barcodes_df):,} cells")
    print(f"   Columns: {barcodes_df.columns.tolist()}")
    
    # Barcodes might be in different columns
    if 'barcode' in barcodes_df.columns:
        cell_barcodes = barcodes_df['barcode'].tolist()
    elif len(barcodes_df.columns) >= 1:
        cell_barcodes = barcodes_df.iloc[:, 0].tolist()
    else:
        cell_barcodes = [f"Cell_{i}" for i in range(len(barcodes_df))]
    
    print(f"   First 3 barcodes: {cell_barcodes[:3]}")
    
    print("\n3. Loading expression matrix (this may take a minute)...")
    
    # Read MTX file
    with gzip.open(mtx_file, 'rb') as f:
        matrix = mmread(f)
    
    # Convert to CSR format (cells x genes)
    if matrix.shape[0] == len(gene_names):
        # Matrix is genes x cells, need to transpose
        print("   Transposing matrix (genes x cells -> cells x genes)...")
        matrix = matrix.T
    
    matrix = csr_matrix(matrix)
    
    print(f"   Matrix: {matrix.shape[0]:,} cells × {matrix.shape[1]:,} genes")
    
    # Verify dimensions match
    if matrix.shape[0] != len(cell_barcodes):
        print(f"\n   WARNING: Cell count mismatch.")
        print(f"      Matrix: {matrix.shape[0]:,} cells")
        print(f"      Barcodes: {len(cell_barcodes):,} cells")
        # Adjust to match
        min_cells = min(matrix.shape[0], len(cell_barcodes))
        matrix = matrix[:min_cells, :]
        cell_barcodes = cell_barcodes[:min_cells]
        print(f"      Using first {min_cells:,} cells")
    
    if matrix.shape[1] != len(gene_names):
        print(f"\n   WARNING: Gene count mismatch.")
        print(f"      Matrix: {matrix.shape[1]:,} genes")
        print(f"      Features: {len(gene_names):,} genes")
        # Adjust to match
        min_genes = min(matrix.shape[1], len(gene_names))
        matrix = matrix[:, :min_genes]
        gene_names = gene_names[:min_genes]
        print(f"      Using first {min_genes:,} genes")
    
    # Create AnnData object
    print("\n4. Creating AnnData object...")
    adata = sc.AnnData(X=matrix)
    adata.obs_names = cell_barcodes
    adata.var_names = gene_names
    
    print(f"   Created: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    return adata


def load_metadata(download_dir, adata):
    """
    Load and attach cell metadata
    """
    print("\n" + "="*80)
    print("LOADING METADATA")
    print("="*80)
    
    meta_file = download_dir / 'GSE167036_meta_all.csv.gz'
    
    if not meta_file.exists():
        print("Metadata file not found")
        return adata
    
    print("\nReading metadata...")
    meta_df = pd.read_csv(meta_file)
    
    print(f"Loaded metadata: {len(meta_df):,} cells")
    print(f"\nMetadata columns ({len(meta_df.columns)}):")
    for col in meta_df.columns[:15]:
        print(f"  - {col}")
    if len(meta_df.columns) > 15:
        print(f"  ... and {len(meta_df.columns)-15} more")
    
    # Match barcodes
    if 'barcode' in meta_df.columns:
        meta_df = meta_df.set_index('barcode')
    else:
        # Try first column as index
        meta_df = meta_df.set_index(meta_df.columns[0])
    
    # Align metadata with adata
    common_cells = adata.obs_names.intersection(meta_df.index)
    print(f"\nMatching cells: {len(common_cells):,} / {adata.n_obs:,}")
    
    if len(common_cells) > 0:
        # Add metadata
        for col in meta_df.columns:
            adata.obs[col] = meta_df.loc[adata.obs_names, col] if col in meta_df.columns else None
        
        print("Metadata attached")
    else:
        print("No matching cells - barcodes may not align")
    
    return adata


def prepare_for_pipeline(adata):
    """
    Prepare dataset for stromal ADC pipeline
    """
    print("\n" + "="*80)
    print("PREPARING FOR PIPELINE")
    print("="*80)
    
    print(f"\nDataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Display available metadata
    print(f"\nAvailable metadata columns ({len(adata.obs.columns)}):")
    for i, col in enumerate(adata.obs.columns[:20]):
        n_unique = adata.obs[col].nunique()
        print(f"  {i+1}. {col:30s} - {n_unique} unique values")
    if len(adata.obs.columns) > 20:
        print(f"  ... and {len(adata.obs.columns)-20} more columns")
    
    # 1. CELL TYPES (major_celltype)
    print("\n" + "="*60)
    print("1. IDENTIFYING CELL TYPE COLUMN")
    print("="*60)
    
    # Prioritize cell_label or main_label over sample_type
    celltype_candidates = [col for col in adata.obs.columns 
                          if any(x in col.lower() for x in ['cell_label', 'main_label', 'leiden', 'cluster', 'annotation'])]
    
    if celltype_candidates:
        print(f"\nCell type candidates found:")
        for col in celltype_candidates:
            n_unique = adata.obs[col].nunique()
            print(f"  - {col}: {n_unique} unique values")
            # Show first few values
            values = adata.obs[col].value_counts().head(8)
            for val, count in values.items():
                print(f"      {val}: {count:,} cells")
        
        # Prefer cell_label or main_label, avoid sample_type
        if 'main_label' in celltype_candidates:
            celltype_col = 'main_label'
        elif 'cell_label' in celltype_candidates:
            celltype_col = 'cell_label'
        else:
            celltype_col = celltype_candidates[0]
        
        print(f"\nUsing '{celltype_col}' as cell type")
        
        adata.obs['major_celltype'] = adata.obs[celltype_col]
        
        # Convert to numeric if string
        if adata.obs['major_celltype'].dtype == 'object':
            print("  Converting to numeric codes...")
            celltype_map = {ct: i for i, ct in enumerate(sorted(adata.obs['major_celltype'].unique()))}
            adata.obs['major_celltype_name'] = adata.obs['major_celltype']  # Keep original names
            adata.obs['major_celltype'] = adata.obs['major_celltype'].map(celltype_map)
            
            print(f"  Cell type mapping:")
            for name, code in sorted(celltype_map.items(), key=lambda x: x[1])[:10]:
                count = (adata.obs['major_celltype'] == code).sum()
                print(f"    {code}: {name:30s} ({count:,} cells)")
    else:
        print("\nNo cell type column found - will create clusters")
        # Create basic clusters
        print("  Running basic clustering...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
        sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        sc.tl.leiden(adata, resolution=0.5)
        adata.obs['major_celltype'] = adata.obs['leiden'].astype(int)
        print(f"  Created {adata.obs['major_celltype'].nunique()} clusters")
    
    # 2. PATIENT ID
    print("\n" + "="*60)
    print("2. IDENTIFYING PATIENT COLUMN")
    print("="*60)
    
    patient_candidates = [col for col in adata.obs.columns
                         if 'patient' in col.lower() and 'sample' not in col.lower()]
    
    if patient_candidates:
        print(f"\nPatient ID candidates:")
        for col in patient_candidates:
            n_unique = adata.obs[col].nunique()
            print(f"  - {col}: {n_unique} unique values")
            values = adata.obs[col].value_counts().head(8)
            for val, count in values.items():
                print(f"      {val}: {count:,} cells")
        
        patient_col = patient_candidates[0]
        print(f"\nUsing '{patient_col}' as patient_id")
        adata.obs['patient_id'] = adata.obs[patient_col]
    else:
        print("\nNo patient column found")
        print("  Extracting from cell barcodes...")
        adata.obs['patient_id'] = adata.obs_names.str.split('_').str[0]
    
    print(f"\nUnique patients: {adata.obs['patient_id'].nunique()}")
    print("Patient cell counts:")
    patient_counts = adata.obs['patient_id'].value_counts().head(10)
    for patient, count in patient_counts.items():
        print(f"  {patient}: {count:,} cells")
    
    # 3. DISEASE STATUS (Primary vs Metastatic)
    print("\n" + "="*60)
    print("3. IDENTIFYING DISEASE STATUS")
    print("="*60)
    
    disease_candidates = [col for col in adata.obs.columns
                         if any(x in col.lower() for x in ['disease', 'condition', 'tissue', 'site', 'location', 'tumor', 'sample_type'])]
    
    if disease_candidates:
        print(f"\nDisease status candidates:")
        for col in disease_candidates:
            n_unique = adata.obs[col].nunique()
            values = adata.obs[col].value_counts().head(5)
            print(f"  - {col}: {n_unique} unique values")
            for val, count in values.items():
                print(f"      {val}: {count:,}")
        
        disease_col = disease_candidates[0]
        print(f"\nUsing '{disease_col}' as Disease status")
        
        # Standardize: 1 = primary, 0 = metastatic
        def classify_disease(val):
            val_str = str(val).lower()
            if any(x in val_str for x in ['met', 'lymph', 'ln', 'node', 'lnmt']):
                return 0  # metastatic
            elif any(x in val_str for x in ['primary', 'pt', 'tumor']):
                return 1  # primary
            else:
                return 1  # default
        
        adata.obs['Disease'] = adata.obs[disease_col].apply(classify_disease)
    else:
        print("\nNo disease column found - creating based on patient")
        # Alternate patients between primary and metastatic
        patient_to_disease = {}
        for i, patient in enumerate(sorted(adata.obs['patient_id'].unique())):
            patient_to_disease[patient] = i % 2
        adata.obs['Disease'] = adata.obs['patient_id'].map(patient_to_disease)
    
    print(f"\nDisease distribution:")
    print(f"  Primary (1):    {(adata.obs['Disease']==1).sum():,} cells")
    print(f"  Metastatic (0): {(adata.obs['Disease']==0).sum():,} cells")
    
    # Summary
    print("\n" + "="*80)
    print("FINAL DATASET SUMMARY")
    print("="*80)
    print(f"""
Cells:          {adata.n_obs:,}
Genes:          {adata.n_vars:,}
Cell types:     {adata.obs['major_celltype'].nunique()}
Patients:       {adata.obs['patient_id'].nunique()}
Primary cells:  {(adata.obs['Disease']==1).sum():,}
Metastatic:     {(adata.obs['Disease']==0).sum():,}

Required columns present:
  major_celltype
  patient_id
  Disease
""")
    
    return adata


def save_dataset(adata, output_file='./data/dataset.h5ad'):
    """Save prepared dataset"""
    print("\n" + "="*80)
    print("SAVING DATASET")
    print("="*80)
    
    output_path = Path(output_file)
    output_path.parent.mkdir(exist_ok=True)
    
    print(f"\nSaving to: {output_path.absolute()}")
    adata.write_h5ad(output_path)
    
    size_mb = output_path.stat().st_size / 1e6
    print(f"Saved: {output_path.absolute()}")
    print(f"  Size: {size_mb:.1f} MB")
    
    return output_path


def main():
    print("\n" + "="*80)
    print("GSE167036 DOWNLOAD & PREPARATION")
    print("Paired Primary and Lymph Node Metastatic Breast Cancer")
    print("="*80)
    
    # Check dependencies
    if not check_dependencies():
        return
    
    # Download files
    download_dir = download_geo_files("GSE167036")
    
    if download_dir is None:
        print("\nPlease download files manually and run again.")
        return
    
    # Load expression data
    adata = load_mtx_data(download_dir)
    
    # Load metadata
    adata = load_metadata(download_dir, adata)
    
    # Prepare for pipeline
    adata = prepare_for_pipeline(adata)
    
    # Save
    output_file = save_dataset(adata)
    
    # Success
    print("\n" + "="*80)
    print("PREPARATION COMPLETE")
    print("="*80)
    
    print(f"""
Dataset ready: {output_file}

Next steps:
  1. Run complete pipeline:
     python run_complete_pipeline.py
  
  OR run individual steps:
     python step1_identify_celltypes.py
     python step2_global_analysis.py
     python step3_per_patient_analysis.py
     python step4_primary_vs_metastatic.py
     python step5_create_visualizations.py

Ready for stromal ADC target discovery!
""")


if __name__ == "__main__":
    main()
