"""
Pipeline steps shared by all datasets.

step2_global   : global + per-subtype stromal targets (surface universe)
step3_samples  : per-sample scoring of the top-N global candidates
step4_stage    : primary vs metastatic (and subtype x stage) aggregation
step6_top15    : top-15 per sample, gene frequencies, summary text

All outputs go to cfg['results'].  Standardized file names (identical in every
dataset):
  celltype_mapping.json, surface_genes.txt
  global_stromal_targets.csv, <subtype>_stromal_targets.csv, subtype_comparison.csv
  per_sample_stromal_targets.csv, sample_statistics.csv, universal_stromal_targets.csv
  individual_samples/<sample>_top15_targets.csv
  primary_stromal_targets.csv, metastatic_stromal_targets.csv, primary_vs_metastatic_comparison.csv
  top15_per_sample.csv, top15_gene_frequencies.csv, top15_summary.txt
"""
import json
import re
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from . import config, scoring, surface, loading

TUMOR_STAGES = ('primary', 'metastatic')


def _safe(name):
    return re.sub(r'[^A-Za-z0-9+_.-]+', '_', str(name))


def archive_old_results(cfg, log=print):
    """Move any pre-existing results into results/_archive_<date>/ so stale and
    fresh outputs can never be mixed."""
    rdir = Path(cfg['results'])
    rdir.mkdir(parents=True, exist_ok=True)
    items = [p for p in rdir.iterdir() if not p.name.startswith('_archive')]
    if not items:
        return None
    dest = rdir / f"_archive_{datetime.now().strftime('%Y-%m-%d')}"
    dest.mkdir(exist_ok=True)
    for p in items:
        target = dest / p.name
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
        shutil.move(str(p), str(target))
    log(f"  archived {len(items)} previous result items -> {dest}")
    return dest


def validate_global(df, log=print):
    top20 = df.head(20)['gene'].tolist()
    found = [g for g in config.KNOWN_STROMAL_SURFACE_MARKERS if g in top20]
    epi = [g for g in config.EPITHELIAL_MARKERS if g in top20]
    log(f"  known surface CAF markers in top 20: {len(found)} {found}")
    if epi:
        log(f"  WARNING epithelial markers in top 20: {epi}")
    else:
        log("  no epithelial markers in top 20")
    return {'known_in_top20': found, 'epithelial_in_top20': epi}


# ---------------------------------------------------------------------------
def step2_global(cfg, adata, log=print, force_surface=False):
    log("\n== STEP 2: global & subtype-specific stromal targets ==")
    rdir = Path(cfg['results'])
    genes, meta = surface.get_surface_genes(adata.var_names, rdir, log, force=force_surface)
    Xs, genes = scoring.subset_matrix(adata, genes)
    s = adata.obs['is_stromal'].values
    t = adata.obs['is_tumor'].values

    raw = scoring.score_genes(Xs, genes, s, t)
    glob = scoring.rank_and_filter(raw)
    glob.to_csv(rdir / 'global_stromal_targets.csv', index=False)
    raw.to_csv(rdir / 'global_stromal_targets_unfiltered.csv', index=False)
    log(f"  global: {len(glob)} genes pass detection >= {config.MIN_STROMAL_DETECTION:.0%}")
    for _, r in glob.head(15).iterrows():
        log(f"    {int(r['rank']):2d}. {r['gene']:10s} excl={r['stromal_exclusivity']:7.3f} "
            f"log2fc={r['log2fc_vs_all']:5.2f} det={r['stromal_detection']:.0%} pen={r['tumor_penalty']:.1f}")
    val = validate_global(glob, log)

    # ---- per-subtype (tumor samples only; Normal/preneoplastic kept as their own groups)
    sub_results = {}
    obs = adata.obs
    for st in sorted(obs['subtype'].unique()):
        if st in ('Unknown', 'nan', ''):
            continue
        u = (obs['subtype'] == st).values
        n_s = int((u & s).sum())
        n_samples = obs.loc[u, 'sample_id'].nunique()
        if n_s < config.MIN_STROMAL_CELLS_PER_SUBTYPE or n_samples < config.MIN_SAMPLES_PER_SUBTYPE:
            log(f"  subtype {st}: skipped ({n_s} stromal cells, {n_samples} samples)")
            continue
        df = scoring.rank_and_filter(scoring.score_genes(Xs, genes, s, t, u))
        df.insert(1, 'subtype', st)
        df.to_csv(rdir / f'{_safe(st)}_stromal_targets.csv', index=False)
        sub_results[st] = df
        log(f"  subtype {st}: {n_samples} samples, {n_s:,} stromal cells -> "
            f"top: {', '.join(df.head(5)['gene'])}")

    if len(sub_results) >= 2:
        comp = compare_subtypes(sub_results)
        comp.to_csv(rdir / 'subtype_comparison.csv', index=False)
        log(f"  subtype_comparison.csv: {len(comp)} genes, classes: "
            + ", ".join(f"{k}={v}" for k, v in comp['target_class'].value_counts().items()))

    with open(rdir / 'step2_summary.json', 'w') as f:
        json.dump({'surface_meta': meta, 'validation': val,
                   'subtypes_analyzed': list(sub_results), 'timestamp': datetime.now().isoformat()},
                  f, indent=2)
    return glob, sub_results


def compare_subtypes(sub_results, threshold=0.5):
    frames = []
    for st, df in sub_results.items():
        frames.append(df[['gene', 'stromal_exclusivity']].assign(subtype=st))
    comb = pd.concat(frames)
    piv = comb.pivot_table(index='gene', columns='subtype', values='stromal_exclusivity',
                           aggfunc='mean').fillna(0)
    subtypes = list(piv.columns)
    piv['mean_exclusivity'] = piv[subtypes].mean(axis=1)
    piv['std_exclusivity'] = piv[subtypes].std(axis=1)
    piv['n_subtypes_high'] = (piv[subtypes] > threshold).sum(axis=1)

    def classify(row):
        high = [c for c in subtypes if row[c] > threshold]
        if len(high) == 0:
            return 'weak'
        if len(high) == 1:
            return f'{high[0]}-specific'
        if len(high) == len(subtypes) or len(high) >= 3:
            return 'pan-subtype'
        return 'multi-subtype'
    piv['target_class'] = piv.apply(classify, axis=1)
    return piv.sort_values('mean_exclusivity', ascending=False).reset_index()


# ---------------------------------------------------------------------------
def step3_samples(cfg, adata, log=print):
    log("\n== STEP 3: per-sample analysis ==")
    rdir = Path(cfg['results'])
    glob = pd.read_csv(rdir / 'global_stromal_targets.csv')
    cand = glob.head(config.N_GLOBAL_CANDIDATES)['gene'].tolist()
    Xs, cand = scoring.subset_matrix(adata, cand)
    log(f"  scoring top-{len(cand)} global candidates in each sample")

    obs = adata.obs
    s = obs['is_stromal'].values
    t = obs['is_tumor'].values
    meta_cols = ['sample_id', 'patient_id', 'stage', 'stage_detail', 'subtype']
    samples = obs[meta_cols].drop_duplicates('sample_id').sort_values('sample_id')

    ind_dir = rdir / 'individual_samples'
    ind_dir.mkdir(exist_ok=True)
    res, stats, skipped = [], [], []
    for _, row in samples.iterrows():
        u = (obs['sample_id'] == row['sample_id']).values
        n_s, n_t, n_u = int((u & s).sum()), int((u & t).sum()), int(u.sum())
        rec = {**row.to_dict(), 'n_total_cells': n_u, 'n_stromal_cells': n_s, 'n_tumor_cells': n_t,
               'pct_stromal': 100 * n_s / n_u, 'pct_tumor': 100 * n_t / n_u}
        if n_s < config.MIN_STROMAL_CELLS_PER_SAMPLE:
            rec['analyzed'] = False
            skipped.append(rec)
            log(f"  {row['sample_id']:35s} SKIP ({n_s} stromal cells)")
            continue
        rec['analyzed'] = True
        df = scoring.score_genes(Xs, cand, s, t, u)
        for c in meta_cols:
            df[c] = row[c]
        res.append(df)
        stats.append(rec)
        top = scoring.rank_and_filter(df).head(config.TOP_N_PER_SAMPLE)
        top[['rank', 'gene', 'stromal_exclusivity', 'log2fc_vs_all', 'stromal_detection',
             'tumor_penalty', 'detection_specificity', 'tumor_mean', 'stromal_mean'] + meta_cols
            ].to_csv(ind_dir / f"{_safe(row['sample_id'])}_top15_targets.csv", index=False)
        log(f"  {row['sample_id']:35s} {n_s:5d} stromal {n_t:5d} tumor  top: "
            + ", ".join(top.head(3)['gene']))

    results_df = pd.concat(res, ignore_index=True)
    stats_df = pd.DataFrame(stats + skipped)
    results_df.to_csv(rdir / 'per_sample_stromal_targets.csv', index=False)
    stats_df.to_csv(rdir / 'sample_statistics.csv', index=False)

    # universal (across analyzed samples)
    filt = results_df[results_df['stromal_detection'] >= config.MIN_STROMAL_DETECTION]
    agg = filt.groupby('gene').agg(exclusivity_mean=('stromal_exclusivity', 'mean'),
                                   exclusivity_std=('stromal_exclusivity', 'std'),
                                   n_samples=('sample_id', 'nunique'),
                                   log2fc_mean=('log2fc_vs_all', 'mean'),
                                   detection_mean=('stromal_detection', 'mean'),
                                   tumor_penalty_mean=('tumor_penalty', 'mean'))
    min_n = max(2, int(0.3 * len(stats)))
    agg = agg[agg['n_samples'] >= min_n].sort_values('exclusivity_mean', ascending=False)
    agg.reset_index().to_csv(rdir / 'universal_stromal_targets.csv', index=False)
    log(f"  {len(stats)} samples analyzed, {len(skipped)} skipped; "
        f"universal (>= {min_n} samples): {', '.join(agg.head(10).index)}")
    return results_df, stats_df


# ---------------------------------------------------------------------------
def _aggregate(df, label, min_samples=2):
    f = df[df['stromal_detection'] >= config.MIN_STROMAL_DETECTION]
    g = f.groupby('gene').agg(exclusivity_mean=('stromal_exclusivity', 'mean'),
                              exclusivity_std=('stromal_exclusivity', 'std'),
                              n_samples=('sample_id', 'nunique'),
                              log2fc_mean=('log2fc_vs_all', 'mean'),
                              detection_mean=('stromal_detection', 'mean'),
                              tumor_penalty_mean=('tumor_penalty', 'mean'),
                              stromal_mean=('stromal_mean', 'mean'),
                              tumor_mean=('tumor_mean', 'mean'))
    g = g[g['n_samples'] >= min_samples].sort_values('exclusivity_mean', ascending=False)
    g['group'] = label
    return g.reset_index()


def _classify_pair(a, b, hi=0.5, lo=0.3):
    if a > hi and b > hi:
        return 'universal'
    if a > hi and b < lo:
        return 'primary_specific'
    if b > hi and a < lo:
        return 'metastatic_specific'
    return 'weak'


def step4_stage(cfg, log=print):
    log("\n== STEP 4: primary vs metastatic ==")
    rdir = Path(cfg['results'])
    res = pd.read_csv(rdir / 'per_sample_stromal_targets.csv')
    stats = pd.read_csv(rdir / 'sample_statistics.csv')
    stats = stats[stats['analyzed']]
    out = {}
    for stage in TUMOR_STAGES:
        ids = stats.loc[stats['stage'] == stage, 'sample_id']
        if len(ids) == 0:
            log(f"  no {stage} samples")
            continue
        agg = _aggregate(res[res['sample_id'].isin(ids)], stage, min_samples=min(2, len(ids)))
        agg.to_csv(rdir / f'{stage}_stromal_targets.csv', index=False)
        out[stage] = agg
        log(f"  {stage}: {len(ids)} samples -> top: {', '.join(agg.head(8)['gene'])}")

    if len(out) == 2:
        p, m = out['primary'], out['metastatic']
        cols = ['gene', 'exclusivity_mean', 'log2fc_mean', 'detection_mean', 'n_samples']
        comp = p[cols].merge(m[cols], on='gene', how='outer', suffixes=('_primary', '_metastatic')).fillna(0)
        comp['exclusivity_diff'] = comp['exclusivity_mean_primary'] - comp['exclusivity_mean_metastatic']
        comp['avg_exclusivity'] = (comp['exclusivity_mean_primary'] + comp['exclusivity_mean_metastatic']) / 2
        comp['target_type'] = [_classify_pair(a, b) for a, b in
                               zip(comp['exclusivity_mean_primary'], comp['exclusivity_mean_metastatic'])]
        comp = comp.sort_values('avg_exclusivity', ascending=False)
        comp.to_csv(rdir / 'primary_vs_metastatic_comparison.csv', index=False)
        log("  classes: " + ", ".join(f"{k}={v}" for k, v in comp['target_type'].value_counts().items()))

        # subtype x stage (only where both arms have >= 2 samples)
        rows = []
        for st in sorted(stats['subtype'].unique()):
            sp = stats[(stats['subtype'] == st) & (stats['stage'] == 'primary')]['sample_id']
            sm = stats[(stats['subtype'] == st) & (stats['stage'] == 'metastatic')]['sample_id']
            if len(sp) < 2 or len(sm) < 2:
                continue
            a = _aggregate(res[res['sample_id'].isin(sp)], 'primary')
            b = _aggregate(res[res['sample_id'].isin(sm)], 'metastatic')
            c = a[cols].merge(b[cols], on='gene', how='outer', suffixes=('_primary', '_metastatic')).fillna(0)
            c.insert(0, 'subtype', st)
            c['target_type'] = [_classify_pair(x, y) for x, y in
                                zip(c['exclusivity_mean_primary'], c['exclusivity_mean_metastatic'])]
            rows.append(c.sort_values('exclusivity_mean_primary', ascending=False))
        if rows:
            pd.concat(rows).to_csv(rdir / 'primary_vs_metastatic_by_subtype.csv', index=False)
            log(f"  subtype x stage comparison written for {len(rows)} subtypes")
    else:
        log("  only one stage present – no comparison file")
    return out


# ---------------------------------------------------------------------------
def step6_top15(cfg, log=print):
    log("\n== STEP 6: top-15 per sample ==")
    rdir = Path(cfg['results'])
    res = pd.read_csv(rdir / 'per_sample_stromal_targets.csv')
    f = res[(res['stromal_detection'] >= config.MIN_STROMAL_DETECTION) & (res['stromal_exclusivity'] > 0)]
    tops = []
    for sid, g in f.groupby('sample_id'):
        g = g.sort_values('stromal_exclusivity', ascending=False).head(config.TOP_N_PER_SAMPLE).copy()
        g['rank'] = np.arange(1, len(g) + 1)
        tops.append(g)
    top = pd.concat(tops, ignore_index=True)
    top.to_csv(rdir / 'top15_per_sample.csv', index=False)

    n = top['sample_id'].nunique()
    freq = top.groupby('gene').agg(n_samples=('sample_id', 'nunique'),
                                   mean_rank=('rank', 'mean'),
                                   avg_exclusivity=('stromal_exclusivity', 'mean')).reset_index()
    freq['pct_samples'] = 100 * freq['n_samples'] / n
    freq['composite_score'] = freq['pct_samples'] / freq['mean_rank']
    freq = freq.sort_values(['n_samples', 'composite_score'], ascending=False)
    freq.to_csv(rdir / 'top15_gene_frequencies.csv', index=False)

    lines = ['=' * 80, f"TOP 15 STROMAL TARGETS PER SAMPLE – {cfg['name']}", '=' * 80, '',
             f"Samples analyzed: {n}", f"Unique genes in top-15 lists: {len(freq)}", '',
             'MOST FREQUENT TARGETS', '-' * 80]
    for i, r in enumerate(freq.head(20).itertuples(), 1):
        lines.append(f"{i:2d}. {r.gene:12s} {r.n_samples:3d}/{n} samples ({r.pct_samples:5.1f}%) | "
                     f"mean rank {r.mean_rank:4.1f} | composite {r.composite_score:5.1f} | "
                     f"avg excl {r.avg_exclusivity:6.3f}")
    lines += ['', '=' * 80, 'PER-SAMPLE TOP 15', '=' * 80]
    for sid, g in top.groupby('sample_id'):
        r0 = g.iloc[0]
        lines.append(f"\n{sid} (patient {r0['patient_id']}, {r0['subtype']}, {r0['stage']}):")
        lines.append('-' * 80)
        for r in g.itertuples():
            lines.append(f"  {int(r.rank):2d}. {r.gene:12s} | excl {r.stromal_exclusivity:7.3f} | "
                         f"log2fc {r.log2fc_vs_all:5.2f} | det {100*r.stromal_detection:4.1f}%")
    (rdir / 'top15_summary.txt').write_text('\n'.join(lines) + '\n')
    log(f"  {n} samples; most frequent: " + ", ".join(freq.head(8)['gene']))
    return top, freq
