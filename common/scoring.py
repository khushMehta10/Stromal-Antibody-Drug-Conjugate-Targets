"""
Stromal-exclusivity score – the single implementation used for every dataset,
for the global, subtype-level and per-sample analyses.

For a gene g and a cell universe U (all cells, or the cells of one subtype /
one sample) with stromal cells S ⊂ U, tumor cells T ⊂ U, non-stromal N = U \\ S:

    stromal_mean      = mean_{S} x_g            (x = log1p CP10k)
    non_stromal_mean  = mean_{N} x_g
    tumor_mean        = mean_{T} x_g
    stromal_detection = frac_{S}(x_g > 0)
    all_detection     = frac_{U}(x_g > 0)

    log2fc            = log2((stromal_mean + p) / (non_stromal_mean + p)),  p = 0.001
    tumor_penalty     = 0.1 if tumor_mean/stromal_mean > 0.5
                        0.5 if                        > 0.2
                        1.0 otherwise            (only if |T| >= 10 and stromal_mean > 0)
    detection_specificity
                      = min(stromal_detection / all_detection, CAP)   if n_detected_in_S >= 10
                      = 0                                              otherwise
    stromal_exclusivity = log2fc * stromal_detection * tumor_penalty * detection_specificity

Change vs. the original per-dataset scripts (2026-11):
  * detection_specificity used a hard floor `all_detection > 0.01 else 0`, which
    creates a cliff (score jumps from 0 to >20 at the 1 % boundary) and let the
    ratio reach ~45x in samples with tiny stromal fractions.  It is now gated on
    an absolute count (>=10 detected stromal cells) and capped at CAP = 10.
  * computed vectorized (all genes at once) instead of one CSR column at a time.
"""
import numpy as np
import pandas as pd
from scipy import sparse

from . import config


def _group_stats(Xs, mask):
    """mean and detection count per column for the rows in mask (CSR input)."""
    n = int(mask.sum())
    if n == 0:
        k = Xs.shape[1]
        return np.zeros(k), np.zeros(k, dtype=int), 0
    sub = Xs[mask]
    mean = np.asarray(sub.mean(axis=0)).ravel()
    det = np.asarray((sub > 0).sum(axis=0)).ravel().astype(int)
    return mean, det, n


def score_genes(Xs, genes, stromal_mask, tumor_mask, universe_mask=None):
    """
    Xs            : sparse (n_cells x n_genes) log-normalized matrix restricted to `genes`
    genes         : list of gene symbols matching Xs columns
    stromal_mask  : boolean (n_cells,)
    tumor_mask    : boolean (n_cells,)
    universe_mask : boolean (n_cells,) – cells considered (default: all)
    Returns a DataFrame with one row per gene.
    """
    if universe_mask is None:
        universe_mask = np.ones(Xs.shape[0], dtype=bool)
    s_mask = stromal_mask & universe_mask
    t_mask = tumor_mask & universe_mask
    n_mask = universe_mask & ~stromal_mask

    s_mean, s_det, n_s = _group_stats(Xs, s_mask)
    t_mean, _, n_t = _group_stats(Xs, t_mask)
    n_mean, _, n_n = _group_stats(Xs, n_mask)
    _, u_det, n_u = _group_stats(Xs, universe_mask)

    p = config.PSEUDOCOUNT
    stromal_detection = s_det / max(n_s, 1)
    all_detection = u_det / max(n_u, 1)
    log2fc = np.log2((s_mean + p) / (n_mean + p))

    penalty = np.ones_like(s_mean)
    if n_t >= config.MIN_TUMOR_CELLS_FOR_PENALTY:
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(s_mean > 0, t_mean / np.where(s_mean > 0, s_mean, 1), 0.0)
        penalty = np.where(ratio > 0.5, 0.1, np.where(ratio > 0.2, 0.5, 1.0))
        penalty = np.where(s_mean > 0, penalty, 1.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        spec = np.where(all_detection > 0, stromal_detection / np.where(all_detection > 0, all_detection, 1), 0.0)
    spec = np.minimum(spec, config.SPECIFICITY_CAP)
    spec = np.where(s_det >= config.MIN_STROMAL_DETECTED_CELLS, spec, 0.0)

    exclusivity = log2fc * stromal_detection * penalty * spec

    return pd.DataFrame({
        'gene': list(genes),
        'stromal_mean': s_mean,
        'stromal_detection': stromal_detection,
        'n_stromal_detected': s_det,
        'tumor_mean': t_mean,
        'non_stromal_mean': n_mean,
        'log2fc_vs_all': log2fc,
        'tumor_penalty': penalty,
        'detection_specificity': spec,
        'stromal_exclusivity': exclusivity,
        'n_stromal_cells': n_s,
        'n_tumor_cells': n_t,
        'n_cells': n_u,
    })


def subset_matrix(adata, genes):
    """CSR matrix of adata.X restricted (and ordered) to `genes` present in adata."""
    present = [g for g in genes if g in adata.var_names]
    idx = adata.var_names.get_indexer(present)
    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    Xs = X[:, idx]
    if not sparse.isspmatrix_csr(Xs):
        Xs = Xs.tocsr()
    return Xs, present


def rank_and_filter(df, min_detection=config.MIN_STROMAL_DETECTION):
    out = df[df['stromal_detection'] >= min_detection].copy()
    out = out.sort_values('stromal_exclusivity', ascending=False).reset_index(drop=True)
    out.insert(0, 'rank', np.arange(1, len(out) + 1))
    return out
