"""Step 5: figures (same set for every dataset; subtype/stage figures only when applicable)."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

from . import config

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 9

SUBTYPE_COLORS = {'ER+': '#e377c2', 'TNBC': '#9467bd', 'HER2+': '#ff7f0e', 'PR+': '#17becf',
                  'Normal': '#2ca02c', 'BRCA1-preneoplastic': '#bcbd22', 'Unknown': '#7f7f7f'}


def _read(rdir, name):
    p = rdir / name
    return pd.read_csv(p) if p.exists() else None


def fig1_top_targets(rdir, fdir, glob):
    top = glob.head(20)
    fig, ax = plt.subplots(figsize=(7, 6))
    colors = ['#d62728' if g in config.KNOWN_STROMAL_SURFACE_MARKERS else '#1f77b4' for g in top['gene']]
    ax.barh(np.arange(len(top)), top['stromal_exclusivity'], color=colors)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top['gene'])
    ax.invert_yaxis()
    ax.set_xlabel('Stromal exclusivity score')
    ax.set_title('Top 20 global stromal-specific surface targets')
    ax.legend(handles=[Patch(color='#d62728', label='Known CAF surface marker'),
                       Patch(color='#1f77b4', label='Other')], loc='lower right')
    fig.tight_layout()
    fig.savefig(fdir / 'fig1_top_targets.pdf')
    plt.close(fig)


def fig2_volcano(rdir, fdir, glob):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(glob['log2fc_vs_all'], glob['stromal_exclusivity'], s=12, c='gray', alpha=0.4)
    top = glob.head(20)
    ax.scatter(top['log2fc_vs_all'], top['stromal_exclusivity'], s=40, c='#d62728', edgecolor='k', lw=0.5)
    for r in top.head(12).itertuples():
        ax.annotate(r.gene, (r.log2fc_vs_all, r.stromal_exclusivity), fontsize=7, xytext=(3, 2),
                    textcoords='offset points')
    ax.set_xlabel('log2 fold change (stromal vs non-stromal)')
    ax.set_ylabel('Stromal exclusivity score')
    ax.set_title('Stromal specificity landscape')
    fig.tight_layout()
    fig.savefig(fdir / 'fig2_volcano.pdf')
    plt.close(fig)


def fig3_sample_heatmap(rdir, fdir, glob, per_sample, stats):
    genes = glob.head(25)['gene'].tolist()
    sub = per_sample[per_sample['gene'].isin(genes)]
    mat = sub.pivot_table(index='gene', columns='sample_id', values='stromal_exclusivity', fill_value=0)
    mat = mat.loc[[g for g in genes if g in mat.index]]
    stats = stats[stats['analyzed']].set_index('sample_id')
    order = stats.sort_values(['subtype', 'stage']).index
    mat = mat[[c for c in order if c in mat.columns]]
    vmax = np.nanpercentile(mat.values, 98) if mat.size else 1
    fig = plt.figure(figsize=(max(8, 0.22 * mat.shape[1] + 3), 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 20], hspace=0.03)
    ax0, ax1 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    cols = [SUBTYPE_COLORS.get(stats.loc[c, 'subtype'], '#7f7f7f') for c in mat.columns]
    ax0.bar(np.arange(len(cols)), 1, color=cols, width=1)
    ax0.set_xlim(-0.5, len(cols) - 0.5)
    ax0.axis('off')
    present = [s for s in SUBTYPE_COLORS if s in set(stats['subtype'])]
    ax0.legend(handles=[Patch(color=SUBTYPE_COLORS[s], label=s) for s in present],
               loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=7)
    sns.heatmap(mat, ax=ax1, cmap='YlOrRd', vmin=0, vmax=vmax,
                cbar_kws={'label': 'Stromal exclusivity'}, linewidths=0.2, linecolor='white')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Gene')
    ax1.set_xticklabels([f"{c} ({stats.loc[c, 'stage'][:1].upper()})" for c in mat.columns],
                        rotation=90, fontsize=6)
    fig.suptitle('Per-sample exclusivity of the top-25 global targets (P=primary, M=metastatic, N=normal)',
                 fontsize=10)
    fig.savefig(fdir / 'fig3_sample_heatmap.pdf', bbox_inches='tight')
    plt.close(fig)


def fig4_primary_vs_met(rdir, fdir, comp):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(comp['exclusivity_mean_primary'], comp['exclusivity_mean_metastatic'], s=14, c='gray', alpha=0.5)
    top = comp.head(25)
    ax.scatter(top['exclusivity_mean_primary'], top['exclusivity_mean_metastatic'], s=40, c='#1f77b4', edgecolor='k', lw=0.5)
    for r in top.itertuples():
        ax.annotate(r.gene, (r.exclusivity_mean_primary, r.exclusivity_mean_metastatic), fontsize=6,
                    xytext=(3, 2), textcoords='offset points')
    m = max(comp['exclusivity_mean_primary'].max(), comp['exclusivity_mean_metastatic'].max())
    ax.plot([0, m], [0, m], 'k--', alpha=0.3)
    ax.set_xlabel('Primary – mean exclusivity')
    ax.set_ylabel('Metastatic – mean exclusivity')
    ax.set_title('Primary vs metastatic stromal targets')
    fig.tight_layout()
    fig.savefig(fdir / 'fig4_primary_vs_metastatic.pdf')
    plt.close(fig)


def fig5_subtypes(rdir, fdir, comp):
    subtypes = [c for c in comp.columns if c not in
                ('gene', 'mean_exclusivity', 'std_exclusivity', 'n_subtypes_high', 'target_class')]
    top = comp.head(30)
    mat = top.set_index('gene')[subtypes]
    fig, ax = plt.subplots(figsize=(1.2 * len(subtypes) + 3, 9))
    sns.heatmap(mat, ax=ax, cmap='YlOrRd', vmin=0, vmax=np.nanpercentile(mat.values, 98),
                cbar_kws={'label': 'Stromal exclusivity'}, linewidths=0.3, linecolor='white')
    ax.set_xlabel('Subtype')
    ax.set_ylabel('Gene')
    ax.set_title('Top-30 targets across subtypes')
    fig.tight_layout()
    fig.savefig(fdir / 'fig5_subtype_comparison.pdf')
    plt.close(fig)


def step5_figures(cfg, log=print):
    log("\n== STEP 5: figures ==")
    rdir = Path(cfg['results'])
    fdir = rdir / 'figures'
    fdir.mkdir(exist_ok=True)
    glob = _read(rdir, 'global_stromal_targets.csv')
    per_sample = _read(rdir, 'per_sample_stromal_targets.csv')
    stats = _read(rdir, 'sample_statistics.csv')
    fig1_top_targets(rdir, fdir, glob)
    fig2_volcano(rdir, fdir, glob)
    fig3_sample_heatmap(rdir, fdir, glob, per_sample, stats)
    made = ['fig1_top_targets', 'fig2_volcano', 'fig3_sample_heatmap']
    comp = _read(rdir, 'primary_vs_metastatic_comparison.csv')
    if comp is not None:
        fig4_primary_vs_met(rdir, fdir, comp)
        made.append('fig4_primary_vs_metastatic')
    sc = _read(rdir, 'subtype_comparison.csv')
    if sc is not None:
        fig5_subtypes(rdir, fdir, sc)
        made.append('fig5_subtype_comparison')
    t = glob.head(30)[['rank', 'gene', 'stromal_exclusivity', 'log2fc_vs_all', 'stromal_detection',
                       'tumor_penalty', 'detection_specificity']].round(3)
    t.to_csv(fdir / 'table1_top_targets.csv', index=False)
    log("  " + ", ".join(made) + ", table1_top_targets")
