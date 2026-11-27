#!/usr/bin/env python3
"""
Set the sample_type labels in the GSE180286 h5ad file
=====================================================
Sets the sample_type column:
- Labels PRIMARY samples as 'primary'
- Labels LYMPH NODE samples as 'metastatic'

Dataset: GSE180286
- 5 patients (P1-P5)
- 15 samples: 5 primary + 10 lymph nodes (2 per patient)
"""

import h5py
import numpy as np
from pathlib import Path

# Primary breast cancer GSM IDs (should be 'primary')
PRIMARY_GSMS = {
    'GSM5457199',  # P1-primary
    'GSM5457202',  # P2-primary
    'GSM5457205',  # P3-primary
    'GSM5457208',  # P4-primary
    'GSM5457211',  # P5-primary
}

# Lymph node GSM IDs (should be 'metastatic')
LYMPH_NODE_GSMS = {
    'GSM5457200',  # P1-lymph node 1
    'GSM5457201',  # P1-lymph node 2
    'GSM5457203',  # P2-lymph node 1
    'GSM5457204',  # P2-lymph node 2
    'GSM5457206',  # P3-lymph node 1
    'GSM5457207',  # P3-lymph node 2
    'GSM5457209',  # P4-lymph node 1
    'GSM5457210',  # P4-lymph node 2
    'GSM5457212',  # P5-lymph node 1
    'GSM5457213',  # P5-lymph node 2
}

def patch_h5ad():
    """Set each cell's sample_type from its GSM ID"""
    
    print("="*80)
    print("SETTING SAMPLE TYPES IN THE GSE180286 H5AD FILE")
    print("="*80)
    
    h5ad_path = Path('./data/dataset.h5ad')
    
    if not h5ad_path.exists():
        print(f"\nERROR: {h5ad_path} not found! Run this script from the GSE180286 folder.")
        return False
    
    print(f"\nFile: {h5ad_path.absolute()}")
    print(f"   Size: {h5ad_path.stat().st_size / 1024**3:.2f} GB")
    
    # Create backup
    backup_path = h5ad_path.parent / f"{h5ad_path.stem}_backup{h5ad_path.suffix}"
    
    if not backup_path.exists():
        print(f"\nCreating backup: {backup_path.name}")
        import shutil
        shutil.copy2(h5ad_path, backup_path)
        print("   Backup created")
    else:
        print(f"\nBackup already exists: {backup_path.name}")
    
    # Open h5ad file for reading and modification
    print("\nLoading metadata...")
    
    with h5py.File(h5ad_path, 'r+') as f:
        obs_group = f['obs']
        
        # Load gsm_id
        print("   Loading gsm_id column...")
        gsm_item = obs_group['gsm_id']
        
        if isinstance(gsm_item, h5py.Group):
            gsm_codes = gsm_item['codes'][:]
            gsm_categories = gsm_item['categories'][:]
            gsm_categories = [x.decode('utf-8') if isinstance(x, bytes) else str(x) 
                             for x in gsm_categories]
            gsm_ids = [gsm_categories[int(i)] for i in gsm_codes]
        else:
            gsm_data = gsm_item[:]
            gsm_ids = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in gsm_data]
        
        print("   Loading sample_type column...")
        sample_type_item = obs_group['sample_type']
        
        if isinstance(sample_type_item, h5py.Group):
            # Categorical
            sample_codes = sample_type_item['codes'][:]
            sample_categories = sample_type_item['categories'][:]
            sample_categories = [x.decode('utf-8') if isinstance(x, bytes) else str(x) 
                                for x in sample_categories]
            
            print(f"\nCurrent sample_type categories: {sample_categories}")
            
            # Ensure both 'primary' and 'metastatic' categories exist
            if 'primary' not in sample_categories:
                print("\nAdding 'primary' category...")
                sample_categories.append('primary')
            
            if 'metastatic' not in sample_categories:
                print("\nAdding 'metastatic' category...")
                sample_categories.append('metastatic')
            
            # Update categories if needed
            sample_type_item['categories'][...] = [x.encode('utf-8') for x in sample_categories]
            
            primary_code = sample_categories.index('primary')
            metastatic_code = sample_categories.index('metastatic')
            
            print(f"\n   Categories updated: {sample_categories}")
            print(f" 'primary' code: {primary_code}")
            print(f" 'metastatic' code: {metastatic_code}")
            
            # Find cells to update
            print("\nFinding cells to update...")
            primary_cells = []
            metastatic_cells = []
            
            for i, gsm_id in enumerate(gsm_ids):
                if gsm_id in PRIMARY_GSMS:
                    primary_cells.append(i)
                elif gsm_id in LYMPH_NODE_GSMS:
                    metastatic_cells.append(i)
            
            print(f"   Found {len(primary_cells):,} cells from PRIMARY samples")
            print(f"   Found {len(metastatic_cells):,} cells from LYMPH NODE samples")
            
            # Update codes
            current_codes = sample_type_item['codes'][:]
            
            if len(primary_cells) > 0:
                print(f"\nUpdating PRIMARY samples to code {primary_code}...")
                for cell_idx in primary_cells:
                    current_codes[cell_idx] = primary_code
                print(f"   Updated {len(primary_cells):,} cells")
            
            if len(metastatic_cells) > 0:
                print(f"\nUpdating LYMPH NODE samples to code {metastatic_code}...")
                for cell_idx in metastatic_cells:
                    current_codes[cell_idx] = metastatic_code
                print(f"   Updated {len(metastatic_cells):,} cells")
            
            # Write back
            sample_type_item['codes'][...] = current_codes
            
            # Show summary by GSM
            print("\nUpdated samples:")
            print("\n   PRIMARY SAMPLES:")
            for gsm in sorted(PRIMARY_GSMS):
                n_cells = sum(1 for g in gsm_ids if g == gsm)
                if n_cells > 0:
                    patient = gsm_to_patient(gsm)
                    print(f"   {gsm} ({patient}): {n_cells:,} cells → primary")
            
            print("\n   METASTATIC SAMPLES (Lymph Nodes):")
            for gsm in sorted(LYMPH_NODE_GSMS):
                n_cells = sum(1 for g in gsm_ids if g == gsm)
                if n_cells > 0:
                    patient = gsm_to_patient(gsm)
                    print(f"   {gsm} ({patient}): {n_cells:,} cells → metastatic")
        
        else:
            print("\nsample_type is not categorical, using direct string update...")
            
            sample_data = sample_type_item[:]
            sample_types = [x.decode('utf-8') if isinstance(x, bytes) else str(x) 
                           for x in sample_data]
            
            primary_updated = 0
            metastatic_updated = 0
            
            for i, gsm_id in enumerate(gsm_ids):
                if gsm_id in PRIMARY_GSMS:
                    sample_types[i] = 'primary'
                    primary_updated += 1
                elif gsm_id in LYMPH_NODE_GSMS:
                    sample_types[i] = 'metastatic'
                    metastatic_updated += 1
            
            # Convert back to bytes and write
            sample_data_new = [x.encode('utf-8') for x in sample_types]
            sample_type_item[...] = sample_data_new
            
            print(f"   Updated {primary_updated:,} cells to 'primary'")
            print(f"   Updated {metastatic_updated:,} cells to 'metastatic'")
    
    print("\n" + "="*80)
    print("SAMPLE TYPES SET")
    print("="*80)
    
    print(f"\nModified: {h5ad_path}")
    print(f"Backup: {backup_path}")
    
    print("\nNext step: python run_GSE180286.py")
    
    return True


def gsm_to_patient(gsm):
    """Convert GSM ID to patient label"""
    gsm_to_patient_map = {
        'GSM5457199': 'P1-primary',
        'GSM5457200': 'P1-LN1',
        'GSM5457201': 'P1-LN2',
        'GSM5457202': 'P2-primary',
        'GSM5457203': 'P2-LN1',
        'GSM5457204': 'P2-LN2',
        'GSM5457205': 'P3-primary',
        'GSM5457206': 'P3-LN1',
        'GSM5457207': 'P3-LN2',
        'GSM5457208': 'P4-primary',
        'GSM5457209': 'P4-LN1',
        'GSM5457210': 'P4-LN2',
        'GSM5457211': 'P5-primary',
        'GSM5457212': 'P5-LN1',
        'GSM5457213': 'P5-LN2',
    }
    return gsm_to_patient_map.get(gsm, 'Unknown')


if __name__ == "__main__":
    success = patch_h5ad()
    
    if success:
        print("\nAll done.")
