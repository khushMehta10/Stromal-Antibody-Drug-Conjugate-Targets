#!/usr/bin/env python
"""
Cross-dataset figures from composite/ (run composite_scoring.py first).

composite/figures/
  fig2_composite_all.pdf            top-25 genes by composite score, all tumor samples
  fig2b_<group>_composite.pdf       same per group (Primary, Metastatic, ERpos, HER2pos, TNBC)
  fig3_presence_absence.pdf         top-25 genes x all samples (present in sample's top-15?), grouped by dataset
  fig4_primary_vs_metastatic.pdf    composite score primary vs metastatic, top genes labeled
  fig5_subtype_panels.pdf           top-10 by composite in ER+ / HER2+ / TNBC
  fig5b_subtype_heatmap.pdf         composite score of the top-20 (All) genes across subtypes
  fig6_rank_vs_frequency.pdf        mean rank vs frequency (All), color = mean exclusivity
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from common import config  # noqa: E402

COMP = ROOT / 'composite'
FIG = COMP / 'figures'
KNOWN = set(config.KNOWN_STROMAL_SURFACE_MARKERS)
BLUE = '#325ea8'
DATASET_COLORS = {'GSE161529': '#1f77b4', 'GSE167036': '#ff7f0e', 'GSE176078': '#2ca02c',
                  'GSE180286': '#d62728', 'GSE305078': '#9467bd', 'GSE306201': '#8c564b'}

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 9
sns.set_style('whitegrid')


def _summary(group):
    return pd.read_csv(COMP / group / 'gene_frequency_summary.csv')


def _color_known(ax):
    for lab in ax.get_yticklabels():
        if lab.get_text() in KNOWN:
            lab.set_color('crimson')
            lab.set_fontweight('bold')


def barplot(group, title, fname, top_n=25, metric='composite_score', xlabel='Composite score (frequency % / mean rank)'):
    df = _summary(group).sort_values(metric, ascending=False).head(top_n)
    n = int(df['n_group'].iloc[0])
    fig, ax = plt.subplots(figsize=(6, 0.28 * top_n + 1.2))
    ax.barh(np.arange(len(df)), df[metric], color=BLUE, zorder=3)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df['gene'])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(f'{title} (n = {n} samples)')
    ax.grid(True, axis='x', alpha=0.3, zorder=0)
    _color_known(ax)
    fig.tight_layout()
    fig.savefig(FIG / fname)
    plt.close(fig)


def presence_absence(top_n=25):
    used = pd.read_csv(COMP / 'samples_used.csv')
    used = used[used['in_All']].copy()
    used['unit'] = used['dataset'] + '|' + used['sample_id'].astype(str)
    tops = []
    for name in config.DATASETS:
        p = config.get(name)['results'] / 'top15_per_sample.csv'
        if p.exists():
            t = pd.read_csv(p)
            t['unit'] = name + '|' + t['sample_id'].astype(str)
            tops.append(t[['unit', 'gene']])
    top = pd.concat(tops)
    top = top[top['unit'].isin(used['unit'])]
    genes = _summary('All').head(top_n)['gene'].tolist()
    used = used.sort_values(['dataset', 'stage', 'sample_id'])
    mat = pd.DataFrame(0, index=genes, columns=used['unit'])
    for u, g in top.groupby('unit')['gene']:
        for gene in g:
            if gene in mat.index:
                mat.loc[gene, u] = 1
    fig = plt.figure(figsize=(14, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 22], hspace=0.02)
    ax0, ax1 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    cols = [DATASET_COLORS[d] for d in used['dataset']]
    ax0.bar(np.arange(len(cols)), 1, color=cols, width=1)
    ax0.set_xlim(-0.5, len(cols) - 0.5)
    ax0.axis('off')
    ax0.legend(handles=[Patch(color=c, label=f"{d} (n={int((used['dataset'] == d).sum())})")
                        for d, c in DATASET_COLORS.items()],
               loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=7)
    sns.heatmap(mat, ax=ax1, cmap=['white', BLUE], cbar=False, linewidths=0.3, linecolor='#dddddd',
                xticklabels=False)
    ax1.set_xlabel(f'Tumor samples (n = {mat.shape[1]}; P = primary, M = metastatic ordered within dataset)')
    ax1.set_ylabel('Gene')
    _color_known(ax1)
    fig.suptitle(f'Presence of the top-{top_n} composite genes in each sample\'s top-15 list', y=0.95)
    fig.savefig(FIG / 'fig3_presence_absence.pdf', bbox_inches='tight')
    plt.close(fig)
    frac = used.groupby('dataset')['stage'].apply(lambda s: (s == 'metastatic').mean())
    return mat


def primary_vs_met(top_n=25):
    p = _summary('Primary').set_index('gene')['composite_score']
    m = _summary('Metastatic').set_index('gene')['composite_score']
    df = pd.concat([p.rename('primary'), m.rename('metastatic')], axis=1).fillna(0)
    df['mean'] = df.mean(axis=1)
    top = df.sort_values('mean', ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df['primary'], df['metastatic'], s=14, c='gray', alpha=0.5)
    ax.scatter(top['primary'], top['metastatic'], s=45, c=BLUE, edgecolor='k', lw=0.5, zorder=3)
    for g, r in top.iterrows():
        ax.annotate(g, (r['primary'], r['metastatic']), fontsize=7, xytext=(3, 2), textcoords='offset points',
                    color='crimson' if g in KNOWN else 'black')
    lim = max(df['primary'].max(), df['metastatic'].max()) * 1.05
    ax.plot([0, lim], [0, lim], 'k--', alpha=0.3)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(f"Primary composite score (n = {int(_summary('Primary')['n_group'].iloc[0])})")
    ax.set_ylabel(f"Metastatic composite score (n = {int(_summary('Metastatic')['n_group'].iloc[0])})")
    ax.set_title('Primary vs metastatic stromal targets')
    fig.tight_layout()
    fig.savefig(FIG / 'fig4_primary_vs_metastatic.pdf')
    plt.close(fig)


def subtype_panels(top_n=10):
    groups = [('ERpos', 'ER+'), ('HER2pos', 'HER2+'), ('TNBC', 'TNBC')]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, (g, label) in zip(axes, groups):
        df = _summary(g).head(top_n)
        ax.barh(np.arange(len(df)), df['composite_score'], color=BLUE, zorder=3)
        ax.set_yticks(np.arange(len(df)))
        ax.set_yticklabels(df['gene'])
        ax.invert_yaxis()
        ax.set_title(f"{label} (n = {int(df['n_group'].iloc[0])})")
        ax.set_xlabel('Composite score')
        ax.grid(True, axis='x', alpha=0.3, zorder=0)
        _color_known(ax)
    fig.suptitle('Top-10 stromal targets per subtype', y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / 'fig5_subtype_panels.pdf', bbox_inches='tight')
    plt.close(fig)

    genes = _summary('All').head(20)['gene'].tolist()
    mat = pd.DataFrame({label: _summary(g).set_index('gene')['composite_score'].reindex(genes).fillna(0)
                        for g, label in groups})
    fig, ax = plt.subplots(figsize=(4.5, 7))
    sns.heatmap(mat, ax=ax, cmap='YlOrRd', annot=True, fmt='.1f', annot_kws={'size': 7},
                cbar_kws={'label': 'Composite score'}, linewidths=0.3, linecolor='white')
    ax.set_ylabel('Gene')
    ax.set_title('Composite score of the top-20 genes by subtype')
    _color_known(ax)
    fig.tight_layout()
    fig.savefig(FIG / 'fig5b_subtype_heatmap.pdf')
    plt.close(fig)


def rank_vs_frequency():
    df = _summary('All')
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc = ax.scatter(df['mean_rank'], df['frequency_pct'], c=df['mean_exclusivity'], cmap='viridis',
                    s=40, edgecolor='k', lw=0.3)
    plt.colorbar(sc, ax=ax, label='Mean stromal exclusivity')
    for r in df.head(20).itertuples():
        ax.annotate(r.gene, (r.mean_rank, r.frequency_pct), fontsize=7, xytext=(3, 2),
                    textcoords='offset points', color='crimson' if r.gene in KNOWN else 'black')
    ax.set_xlabel('Mean within-sample rank (lower is better)')
    ax.set_ylabel('Frequency in top-15 lists (%)')
    ax.set_title(f"All tumor samples (n = {int(df['n_group'].iloc[0])})")
    fig.tight_layout()
    fig.savefig(FIG / 'fig6_rank_vs_frequency.pdf')
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    barplot('All', 'Top-25 stromal targets – all tumor samples', 'fig2_composite_all.pdf')
    barplot('All', 'Top-25 by frequency – all tumor samples', 'fig2_frequency_all.pdf',
            metric='frequency_pct', xlabel='Frequency in top-15 lists (%)')
    for g, label in [('Primary', 'Primary'), ('Metastatic', 'Metastatic'), ('ERpos', 'ER+'),
                     ('HER2pos', 'HER2+'), ('TNBC', 'TNBC')]:
        barplot(g, f'Top-25 stromal targets – {label}', f'fig2b_{g}_composite.pdf')
    presence_absence()
    primary_vs_met()
    subtype_panels()
    rank_vs_frequency()
    print('figures written to', FIG)
    for p in sorted(FIG.iterdir()):
        print('  ', p.name)


if __name__ == '__main__':
    main()
