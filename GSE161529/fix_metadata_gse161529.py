#!/usr/bin/env python
"""
Rewrite the sample metadata in GSE161529/data/dataset.h5ad from the GEO series
matrix instead of parsing it out of file names.

Bugs in the original file-name parser (download_preprocess_gse161529.py):
  * `'ER' in subtype_code` was tested before `== 'HER2'`, and "HER2" contains
    "ER", so all six HER2+ tumors were labeled ER+.
  * for TN-B1-<patient> names the second token ("B1") was taken as the patient,
    merging four different BRCA1-mutant TNBC patients into one patient "B1".
  * B1-<patient> (BRCA1 pre-neoplastic tissue) was labeled TNBC.
  * normal tissue was labeled stage 'primary'.

New obs columns (old ones are kept with an `orig_` prefix for traceability):
  patient_id      GEO "patient" field (e.g. 0064)
  cancer_type     GEO "cancer type" verbatim
  cancer_subtype  ER+ | PR+ | HER2+ | TNBC | Normal | BRCA1-preneoplastic
  brca1_carrier   True for "BRCA1 pre-neoplastic" and "Triple negative BRCA1 tumour"
  cell_population GEO "cell population": Total | Epithelial | Involved lymph node
  stage           primary | metastatic | normal | preneoplastic
  stage_detail    e.g. "primary tumor (Total)", "involved lymph node", "normal (Epithelial-sorted)"
  sample_type     == stage  (kept for backwards compatibility)

Only the obs group of the h5ad is rewritten (X is untouched).  The previous obs
table is saved to data/obs_before_metadata_fix.csv.
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
H5AD = HERE / 'data' / 'dataset.h5ad'
SERIES = HERE / 'GSE161529_series_matrix.txt'

try:
    from anndata.io import read_elem, write_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem, write_elem


def parse_series_matrix(path):
    rows = {}
    for line in open(path):
        if not line.startswith('!Sample_'):
            continue
        key, *vals = line.rstrip('\n').split('\t')
        vals = [v.strip('"') for v in vals]
        if key == '!Sample_characteristics_ch1':
            k = vals[0].split(':')[0]
            rows[k] = [v.split(': ', 1)[1] if ': ' in v else '' for v in vals]
        elif key in ('!Sample_title', '!Sample_geo_accession'):
            rows[key.replace('!Sample_', '')] = vals
    df = pd.DataFrame(rows).set_index('geo_accession')
    return df


SUBTYPE = {
    'Normal': 'Normal',
    'BRCA1 pre-neoplastic': 'BRCA1-preneoplastic',
    'Triple negative tumour': 'TNBC',
    'Triple negative BRCA1 tumour': 'TNBC',
    'HER2+ tumour': 'HER2+',
    'ER+ tumour': 'ER+',
    'PR+ tumour': 'PR+',
}


def derive(meta_row):
    ct = meta_row['cancer type']
    pop = meta_row['cell population']
    subtype = SUBTYPE[ct]
    brca1 = 'BRCA1' in ct
    if ct == 'Normal':
        stage = 'normal'
        detail = f'normal ({pop}-sorted)' if pop == 'Epithelial' else 'normal (Total)'
    elif ct == 'BRCA1 pre-neoplastic':
        stage = 'preneoplastic'
        detail = 'BRCA1 pre-neoplastic (Total)'
    elif pop == 'Involved lymph node':
        stage = 'metastatic'
        detail = 'involved lymph node'
    else:
        stage = 'primary'
        detail = f'primary tumor ({pop})'
    return dict(patient_id=meta_row['patient'], cancer_type=ct, cancer_subtype=subtype,
                brca1_carrier=brca1, cell_population=pop, stage=stage, stage_detail=detail,
                sample_type=stage)


def main(dry_run=False):
    meta = parse_series_matrix(SERIES)
    print(f"series matrix: {len(meta)} samples")

    with h5py.File(H5AD, 'r') as f:
        obs = read_elem(f['obs'])
    print(f"h5ad obs: {obs.shape[0]:,} cells, columns {list(obs.columns)}")

    backup = HERE / 'data' / 'obs_before_metadata_fix.csv'
    if not backup.exists():
        obs.to_csv(backup)
        print(f"saved previous obs -> {backup}")

    gsm = obs['gsm_id'].astype(str)
    missing = sorted(set(gsm) - set(meta.index))
    if missing:
        sys.exit(f"GSM ids in h5ad but not in series matrix: {missing}")

    per_gsm = {g: derive(meta.loc[g]) for g in gsm.unique()}
    new = pd.DataFrame([per_gsm[g] for g in gsm], index=obs.index)

    for col in ['patient_id', 'cancer_subtype', 'sample_type']:
        obs[f'orig_{col}'] = obs[col].astype(str).values
    for col in new.columns:
        obs[col] = new[col].values
    for col in ['patient_id', 'cancer_type', 'cancer_subtype', 'cell_population',
                'stage', 'stage_detail', 'sample_type', 'orig_patient_id',
                'orig_cancer_subtype', 'orig_sample_type']:
        obs[col] = obs[col].astype('category')
    obs['brca1_carrier'] = obs['brca1_carrier'].astype(bool)

    summary = (obs.drop_duplicates('sample_id')
                  .groupby(['cancer_subtype', 'stage'], observed=True)['sample_id'].count())
    print("\nsamples per subtype x stage:\n" + summary.to_string())
    changed = (obs['orig_cancer_subtype'].astype(str) != obs['cancer_subtype'].astype(str))
    print(f"\ncells whose subtype label changed: {changed.sum():,} "
          f"({obs.loc[changed, 'sample_id'].nunique()} samples)")
    changed_p = (obs['orig_patient_id'].astype(str) != obs['patient_id'].astype(str))
    print(f"cells whose patient id changed: {changed_p.sum():,} "
          f"({obs.loc[changed_p, 'sample_id'].nunique()} samples)")

    if dry_run:
        print("\n--dry-run: h5ad not modified")
        return
    with h5py.File(H5AD, 'r+') as f:
        del f['obs']
        write_elem(f, 'obs', obs)
    print(f"\nwrote new obs into {H5AD}")


if __name__ == '__main__':
    main(dry_run='--dry-run' in sys.argv)
