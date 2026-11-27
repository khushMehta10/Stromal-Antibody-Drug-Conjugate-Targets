#!/usr/bin/env python
"""
Cross-dataset composite scoring of stromal ADC targets.

Pools every dataset's per-sample top-15 lists and computes, per gene and per
group (All / Primary / Metastatic / ER+ / HER2+ / TNBC):

    n_samples       number of samples whose top-15 contains the gene
    frequency_pct   100 * n_samples / n_samples_in_group
    mean_rank       mean within-sample rank (1 = best)
    composite_score frequency_pct / mean_rank          <- the ranking score

This reproduces the formula of ORIG_GFA/all_folder/gene_frequency_analysis.py but
builds the groups automatically from each dataset's sample_statistics.csv
instead of hand-copied folders.

Group membership
  All         : tumor samples only (stage primary or metastatic); normal and
                BRCA1 pre-neoplastic samples are excluded
  Primary / Metastatic : by stage
  ER+ / HER2+ / TNBC   : tumor samples by subtype (PR+ and Unknown excluded)

Outputs (composite/):
  <group>/gene_frequency_summary.csv       full table for the group
  <group>/top25_composite.txt              gene list (one per line) for the bulk validation scripts
  composite_summary.csv                    composite score of every gene in every group (wide)
  sample_counts.csv                        samples per dataset x group (the n reported in the text)
  samples_used.csv                         every sample with dataset / stage / subtype / group flags

    python composite_scoring.py                 # current results of all six datasets
    python composite_scoring.py --legacy ORIG_GFA   # recompute from the original hand-made folders
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from common import config  # noqa: E402

OUT = ROOT / 'composite'
TOP_N_LIST = 25
SUBTYPES = ['ER+', 'HER2+', 'TNBC']


def summarize(top_df, unit_col='sample_id'):
    n = top_df[unit_col].nunique()
    g = top_df.groupby('gene').agg(n_samples=(unit_col, 'nunique'),
                                   mean_rank=('rank', 'mean'),
                                   median_rank=('rank', 'median'),
                                   mean_exclusivity=('stromal_exclusivity', 'mean'),
                                   mean_log2fc=('log2fc_vs_all', 'mean')).reset_index()
    g['n_group'] = n
    g['frequency_pct'] = 100 * g['n_samples'] / n
    g['composite_score'] = g['frequency_pct'] / g['mean_rank']
    g['n_datasets'] = top_df.groupby('gene')['dataset'].nunique().reindex(g['gene']).values \
        if 'dataset' in top_df.columns else np.nan
    return g.sort_values(['composite_score', 'frequency_pct'], ascending=False).reset_index(drop=True)


def load_current():
    frames, stats = [], []
    for name in config.DATASETS:
        rdir = config.get(name)['results']
        t = rdir / 'top15_per_sample.csv'
        s = rdir / 'sample_statistics.csv'
        if not t.exists():
            print(f"  {name}: no top15_per_sample.csv – skipped")
            continue
        top = pd.read_csv(t)
        top['dataset'] = name
        top['unit'] = name + '|' + top['sample_id'].astype(str)
        frames.append(top)
        st = pd.read_csv(s)
        st = st[st['analyzed']].copy()
        st['dataset'] = name
        st['unit'] = name + '|' + st['sample_id'].astype(str)
        stats.append(st)
    return pd.concat(frames, ignore_index=True), pd.concat(stats, ignore_index=True)


def groups_from_stats(stats):
    tum = stats['stage'].isin(['primary', 'metastatic'])
    g = {'All': stats[tum],
         'Primary': stats[stats['stage'] == 'primary'],
         'Metastatic': stats[stats['stage'] == 'metastatic']}
    for st in SUBTYPES:
        g[st] = stats[tum & (stats['subtype'] == st)]
    return g


def run_current():
    top, stats = load_current()
    OUT.mkdir(exist_ok=True)
    groups = groups_from_stats(stats)
    wide = None
    counts = []
    flags = stats[['dataset', 'sample_id', 'patient_id', 'stage', 'stage_detail', 'subtype',
                   'n_stromal_cells', 'n_tumor_cells']].copy()
    for gname, gstats in groups.items():
        units = set(gstats['unit'])
        flags[f'in_{gname}'] = stats['unit'].isin(units).values
        sub = top[top['unit'].isin(units)]
        if sub.empty:
            continue
        summ = summarize(sub, unit_col='unit')
        d = OUT / gname.replace('+', 'pos')
        d.mkdir(exist_ok=True)
        summ.to_csv(d / 'gene_frequency_summary.csv', index=False)
        (d / f'top{TOP_N_LIST}_composite.txt').write_text('\n'.join(summ.head(TOP_N_LIST)['gene']) + '\n')
        col = summ.set_index('gene')['composite_score'].rename(gname)
        wide = col.to_frame() if wide is None else wide.join(col, how='outer')
        per_ds = gstats.groupby('dataset')['unit'].nunique()
        counts.append(per_ds.rename(gname))
        print(f"\n{gname}: {len(units)} samples "
              + "(" + ", ".join(f"{k}={v}" for k, v in per_ds.items()) + ")")
        for r in summ.head(10).itertuples():
            print(f"  {r.gene:10s} composite={r.composite_score:5.1f}  freq={r.frequency_pct:5.1f}%  "
                  f"mean_rank={r.mean_rank:4.1f}  n={r.n_samples}/{r.n_group}  datasets={r.n_datasets}")
    wide = wide.fillna(0).sort_values('All', ascending=False)
    wide.to_csv(OUT / 'composite_summary.csv')
    cnt = pd.concat(counts, axis=1).fillna(0).astype(int)
    cnt.loc['TOTAL'] = cnt.sum()
    cnt.to_csv(OUT / 'sample_counts.csv')
    flags.to_csv(OUT / 'samples_used.csv', index=False)
    print("\nSamples per dataset x group:\n" + cnt.to_string())
    print(f"\nwritten to {OUT}/")


def run_legacy(base):
    """Recompute the original hand-made analysis (folders of per-sample CSVs)."""
    base = Path(base)
    folders = {'All': base / 'all_folder' / 'All',
               'Primary': base / 'all_stage_folder' / 'Primary',
               'Metastatic': base / 'all_stage_folder' / 'Metastatic',
               'ER+': base / 'all_subtype_folder' / 'ER+',
               'HER2+': base / 'all_subtype_folder' / 'HER2+',
               'TNBC': base / 'all_subtype_folder' / 'TNBC'}
    for gname, d in folders.items():
        files = sorted(d.glob('*.csv'))
        frames = []
        for f in files:
            df = pd.read_csv(f)
            df['sample_id'] = f.stem
            frames.append(df)
        top = pd.concat(frames)
        summ = summarize(top)
        print(f"\n[legacy] {gname}: {len(files)} samples")
        for r in summ.head(8).itertuples():
            print(f"  {r.gene:10s} composite={r.composite_score:5.1f}  freq={r.frequency_pct:5.1f}%  mean_rank={r.mean_rank:4.1f}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--legacy', help='path to ORIG_GFA to recompute the original analysis')
    a = ap.parse_args()
    if a.legacy:
        run_legacy(a.legacy)
    else:
        run_current()
