"""
Load a dataset's h5ad and standardize it for the engine.

Guarantees after load_dataset():
  * adata.X is log1p(CP10k)-normalized (three of the six h5ads stored raw UMI
    counts in X; the old per-dataset scripts scored those on raw counts).
  * adata.obs has: sample_id, patient_id, stage, stage_detail, subtype,
    is_stromal, is_tumor
"""
import json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from . import config


def _looks_like_raw_counts(X, n_check=200_000):
    if sparse.issparse(X):
        vals = X.data[:n_check]
    else:
        vals = np.asarray(X).ravel()[:n_check]
    vals = vals[vals != 0]
    if vals.size == 0:
        return False
    is_int = np.allclose(vals, np.round(vals))
    return bool(is_int and vals.max() > 20)


def ensure_lognorm(adata, log):
    if _looks_like_raw_counts(adata.X):
        log(f"  X holds raw UMI counts (max={float(adata.X.max()):.0f}) -> normalize_total(1e4) + log1p")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns['normalized_by_engine'] = True
    else:
        log(f"  X already log-normalized (max={float(adata.X.max()):.2f})")
        adata.uns['normalized_by_engine'] = False
    return adata


def _read_h5ad(path, log):
    adata = sc.read_h5ad(path)
    # GSE306201 (Seurat->h5ad conversion) stores categorical obs as int codes
    # with the category strings in obs['__categories'] (an h5 group that anndata
    # does not map to a column). Decode them here.
    import h5py
    with h5py.File(path, 'r') as f:
        if 'obs' in f and '__categories' in f['obs']:
            grp = f['obs']['__categories']
            for col in grp.keys():
                if col in adata.obs.columns:
                    cats = [c.decode() if isinstance(c, bytes) else str(c) for c in grp[col][:]]
                    codes = adata.obs[col].astype(int).values
                    adata.obs[col] = pd.Categorical.from_codes(codes, cats)
            log(f"  decoded {len(grp.keys())} coded obs columns via obs/__categories")
    return adata


def load_dataset(cfg, log=print):
    log(f"Loading {cfg['h5ad']}")
    adata = _read_h5ad(cfg['h5ad'], log)
    log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    ensure_lognorm(adata, log)

    obs = adata.obs
    ct = obs[cfg['celltype_col']].astype(str)
    obs['is_stromal'] = ct.isin([str(x) for x in cfg['stromal_labels']]).values
    obs['is_tumor'] = ct.isin([str(x) for x in cfg['tumor_labels']]).values

    obs['sample_id'] = obs[cfg['sample_col']].astype(str).values
    obs['patient_id'] = obs[cfg['patient_col']].astype(str).values

    # ---- stage -----------------------------------------------------------
    if 'stage_col' in cfg:
        obs['stage'] = obs[cfg['stage_col']].astype(str).values
    elif 'stage_map' in cfg:
        col, m = cfg['stage_map']
        obs['stage'] = obs[col].astype(str).map(m).fillna('unknown').values
    elif 'stage_const' in cfg:
        obs['stage'] = cfg['stage_const']
    else:
        raise ValueError('no stage definition')

    if 'stage_detail_col' in cfg:
        obs['stage_detail'] = obs[cfg['stage_detail_col']].astype(str).values
    elif 'stage_detail_map' in cfg:
        col, m = cfg['stage_detail_map']
        obs['stage_detail'] = obs[col].astype(str).map(m).fillna('').values
    elif 'stage_detail' in obs.columns:
        obs['stage_detail'] = obs['stage_detail'].astype(str).values
    else:
        obs['stage_detail'] = obs['stage'].values

    # ---- subtype ---------------------------------------------------------
    if 'subtype_col' in cfg:
        obs['subtype'] = obs[cfg['subtype_col']].astype(str).values
    elif 'subtype_fn' in cfg:
        obs['subtype'] = cfg['subtype_fn'](obs).astype(str).values
    elif 'subtype_const' in cfg:
        obs['subtype'] = cfg['subtype_const']
    else:
        obs['subtype'] = 'Unknown'

    n_s, n_t = int(obs['is_stromal'].sum()), int(obs['is_tumor'].sum())
    log(f"  stromal ({cfg['stromal_labels']}): {n_s:,} ({100*n_s/adata.n_obs:.1f}%)")
    log(f"  tumor   ({cfg['tumor_labels']}): {n_t:,} ({100*n_t/adata.n_obs:.1f}%)")
    if n_s == 0 or n_t == 0:
        raise ValueError('stromal or tumor label set matched zero cells – check config')
    log("  stage:   " + ", ".join(f"{k}={v}" for k, v in obs['stage'].value_counts().items()))
    log("  subtype: " + ", ".join(f"{k}={v}" for k, v in obs['subtype'].value_counts().items()))
    return adata


def write_celltype_mapping(cfg, adata):
    """Record the definitions actually used (replaces the old step-1 output)."""
    obs = adata.obs
    mapping = {
        'dataset': cfg['name'],
        'celltype_column': cfg['celltype_col'],
        'stromal_celltypes': [str(x) for x in cfg['stromal_labels']],
        'tumor_celltypes': [str(x) for x in cfg['tumor_labels']],
        'n_cells': int(adata.n_obs),
        'n_stromal': int(obs['is_stromal'].sum()),
        'n_tumor': int(obs['is_tumor'].sum()),
        'pct_stromal': float(100 * obs['is_stromal'].mean()),
        'pct_tumor': float(100 * obs['is_tumor'].mean()),
        'normalized_by_engine': bool(adata.uns.get('normalized_by_engine', False)),
        'stromal_definition': 'fibroblast/CAF lineage only (no endothelial)',
        'method': 'annotation-based (see common/config.py)',
    }
    out = cfg['results'] / 'celltype_mapping.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w') as f:
        json.dump(mapping, f, indent=2)
    return mapping
