#!/usr/bin/env python
"""
Produce the figures with the ORIGINAL figure code (the step5 plotting scripts in
original_plots/ and the ORIG_GFA gene-frequency scripts), fed with the new pipeline results.

Per dataset
  * builds the input files the original step5 expects (legacy names/columns)
    in  <dataset>/results/_legacy_view/results/
  * runs original_plots/<dataset>/step5_*.py there
  * moves the PDFs into <dataset>/results/figures/   (engine figures -> figures_engine/)

Composite
  * builds composite/original_style/{all_folder,all_stage_folder,all_subtype_folder}
    with one per-sample top-15 CSV per sample, exactly like the hand-made ORIG_GFA
    folders, and runs the three ORIG_GFA gene_frequency_analysis scripts there.

Edits made to the copies of the original scripts (the originals are not changed):
  * step5: the commented-out fig3 heatmap call (GSE161529, GSE180286) is re-enabled
  * GFA: units are labeled samples, not patients, and stage panels say "Per-Stage"
  * GFA: red gene labels mark known ADC targets (KNOWN_ADC_TARGETS below) instead of the
    scripts' hand-picked list of classic fibroblast markers
  * GFA: the presence/absence heatmap no longer prints sample IDs under its columns (100+
    overlapping labels were unreadable) and labels the axis "Tumor samples (n = N)"
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from common import config  # noqa: E402

PLOT_SCRIPTS = ROOT / 'original_plots'
PY = sys.executable

STEP5 = {
    'GSE161529': 'step5_subtype_viz.py',
    'GSE167036': 'step5_create_visualizations.py',
    'GSE176078': 'step5_subtypes_viz.py',
    'GSE180286': 'step5_create_visualizations.py',
    'GSE305078': 'step5_create_visualizations.py',
    'GSE306201': 'step5_create_visualizations.py',
}


# Red gene labels in the composite figures = known ADC targets: genes whose protein is
# targeted by an antibody-drug conjugate that has entered clinical trials (checked
# 2026-11-16; each entry below gives its source). Only the genes shown in the figures
# were checked; ANTXR1, CD99 and others with preclinical-only ADCs are deliberately
# not included. Add a gene here only with a source.
KNOWN_ADC_TARGETS = {
    'EGFR': 'depatuxizumab mafodotin, phase III (Lassman et al., 2023, doi:10.1093/neuonc/noac173)',
    'FAP': 'OMTX705, phase 1 (NCT05547321)',
    'GPC3': 'ZW251, phase 1 (NCT07164313)',
    'LRRC15': 'ABBV-085 / samrotamab vedotin, phase 1 (Demetri et al., 2021, doi:10.1158/1078-0432.CCR-20-4513)',
    'MCAM': 'AMT-253, anti-MUC18, phase 1/2 (NCT06209580; Shi et al., 2023, doi:10.1158/0008-5472.CAN-23-1356)',
}


def _copy_if(src, dst):
    if src.exists():
        shutil.copy(src, dst)


def build_legacy_inputs(name):
    cfg = config.get(name)
    rdir = cfg['results']
    view = rdir / '_legacy_view' / 'results'
    if view.parent.exists():
        shutil.rmtree(view.parent)
    view.mkdir(parents=True)

    for f in ['global_stromal_targets.csv', 'universal_stromal_targets.csv',
              'primary_stromal_targets.csv', 'metastatic_stromal_targets.csv',
              'primary_vs_metastatic_comparison.csv', 'primary_vs_metastatic_by_subtype.csv']:
        _copy_if(rdir / f, view / f)
    for f in rdir.glob('*_stromal_targets.csv'):
        _copy_if(f, view / f.name)

    per = pd.read_csv(rdir / 'per_sample_stromal_targets.csv')
    stats = pd.read_csv(rdir / 'sample_statistics.csv')
    stats = stats[stats['analyzed']].copy()

    # legacy column aliases
    per['cancer_subtype'] = per['subtype']
    stats['cancer_subtype'] = stats['subtype']
    stats['disease_status'] = stats['stage']
    if name == 'GSE176078':
        # the old GSE176078 step3 stored the SUBTYPE in sample_type
        per['sample_type'] = per['subtype']
        stats['sample_type'] = stats['subtype']
    else:
        per['sample_type'] = per['stage']
        stats['sample_type'] = stats['stage']
    stats['gsm_id'] = stats['sample_id']

    # the old "per_patient" tables were really per-sample tables; keep that
    per_pat = per.copy()
    per_pat['patient_id'] = per_pat['sample_id']
    stats_pat = stats.copy()
    stats_pat['patient_id'] = stats_pat['sample_id']

    per.to_csv(view / 'per_sample_stromal_targets.csv', index=False)
    stats.to_csv(view / 'sample_statistics.csv', index=False)
    per_pat.to_csv(view / 'per_patient_stromal_targets.csv', index=False)
    stats_pat.to_csv(view / 'patient_statistics.csv', index=False)

    if name == 'GSE176078':
        # old step4 format: ER_/HER2_/TNBC_stromal_targets.csv + subtype_comparison.csv
        frames = {}
        for st, tag in [('ER+', 'ER'), ('HER2+', 'HER2'), ('TNBC', 'TNBC')]:
            f = rdir / f'{st}_stromal_targets.csv'
            if not f.exists():
                continue
            df = pd.read_csv(f)
            agg = df[['gene', 'stromal_exclusivity', 'log2fc_vs_all', 'stromal_detection']].copy()
            agg.columns = ['gene', 'exclusivity_mean', 'log2fc_mean', 'detection_mean']
            agg['n_patients'] = int(stats.loc[stats['subtype'] == st, 'sample_id'].nunique())
            agg['subtype'] = st
            agg.to_csv(view / f'{tag}_stromal_targets.csv', index=False)
            frames[tag] = agg
        if len(frames) == 3:
            comp = frames['ER'][['gene', 'exclusivity_mean', 'log2fc_mean', 'detection_mean', 'n_patients']].merge(
                frames['HER2'][['gene', 'exclusivity_mean', 'log2fc_mean', 'detection_mean', 'n_patients']],
                on='gene', suffixes=('_ER', '_HER2'), how='outer').merge(
                frames['TNBC'][['gene', 'exclusivity_mean', 'log2fc_mean', 'detection_mean', 'n_patients']],
                on='gene', how='outer').rename(columns={
                    'exclusivity_mean': 'exclusivity_mean_TNBC', 'log2fc_mean': 'log2fc_mean_TNBC',
                    'detection_mean': 'detection_mean_TNBC', 'n_patients': 'n_patients_TNBC'}).fillna(0)
            e = comp[['exclusivity_mean_ER', 'exclusivity_mean_HER2', 'exclusivity_mean_TNBC']]
            comp['avg_exclusivity'] = e.mean(axis=1)
            comp['exclusivity_cv'] = e.std(axis=1, ddof=0) / (e.mean(axis=1) + 0.001)

            def classify(r):
                er, h, t, cv = r['exclusivity_mean_ER'], r['exclusivity_mean_HER2'], r['exclusivity_mean_TNBC'], r['exclusivity_cv']
                if er > 30 and h > 30 and t > 30 and cv < 0.5:
                    return 'universal'
                if er > 30 and h < 20 and t < 20:
                    return 'ER_specific'
                if h > 30 and er < 20 and t < 20:
                    return 'HER2_specific'
                if t > 30 and er < 20 and h < 20:
                    return 'TNBC_specific'
                if max(er, h, t) > 20:
                    return 'moderate'
                return 'weak'
            comp['target_type'] = comp.apply(classify, axis=1)
            comp.sort_values('avg_exclusivity', ascending=False).to_csv(view / 'subtype_comparison.csv', index=False)
    else:
        _copy_if(rdir / 'subtype_comparison.csv', view / 'subtype_comparison.csv')
    return view.parent


def run_dataset(name):
    cfg = config.get(name)
    rdir = cfg['results']
    work = build_legacy_inputs(name)
    src = PLOT_SCRIPTS / name / STEP5[name]
    script = work / STEP5[name]
    text = src.read_text()
    # GSE161529 and GSE180286 had the fig3 heatmap call commented out
    text = text.replace('    # create_patient_heatmap(data, output_dir)', '    create_patient_heatmap(data, output_dir)')
    text = text.replace('        print(matrix)\n', '')
    script.write_text(text)
    print(f"\n=== {name}: running original {STEP5[name]}")
    r = subprocess.run([PY, script.name], cwd=work, capture_output=True, text=True)
    tail = '\n'.join(r.stdout.strip().splitlines()[-25:])
    print(tail)
    if r.returncode != 0:
        print(r.stderr[-3000:])
        raise SystemExit(f"{name} step5 failed")

    fig_src = work / 'results' / 'figures'
    fig_dst = rdir / 'figures'
    eng = rdir / 'figures_engine'
    if fig_dst.exists() and not eng.exists():
        shutil.move(str(fig_dst), str(eng))
    elif fig_dst.exists():
        shutil.rmtree(fig_dst)
    shutil.move(str(fig_src), str(fig_dst))
    print(f"  -> {fig_dst}: {', '.join(sorted(p.name for p in fig_dst.iterdir()))}")


# ---------------------------------------------------------------------------
def build_composite_folders():
    base = ROOT / 'composite' / 'original_style'
    if base.exists():
        shutil.rmtree(base)
    layout = {'all_folder': ['All'], 'all_stage_folder': ['Primary', 'Metastatic'],
              'all_subtype_folder': ['ER+', 'HER2+', 'TNBC']}
    for folder, groups in layout.items():
        for g in groups:
            (base / folder / g).mkdir(parents=True)
    counts = {g: 0 for gs in layout.values() for g in gs}
    # Per-sample lists come from top15_per_sample.csv (step 6), which keeps only
    # genes with stromal_exclusivity > 0. The step-3 files in individual_samples/
    # pad a sample's list with zero-score genes when fewer than 15 genes score
    # above 0, which made these figures disagree slightly with composite_scoring.py.
    cols = ['rank', 'gene', 'stromal_exclusivity', 'log2fc_vs_all', 'stromal_detection', 'tumor_penalty',
            'detection_specificity', 'tumor_mean', 'stromal_mean', 'sample_id', 'patient_id', 'stage',
            'stage_detail', 'subtype']
    for name in config.DATASETS:
        rdir = config.get(name)['results']
        stats = pd.read_csv(rdir / 'sample_statistics.csv')
        stats = stats[stats['analyzed']]
        top = pd.read_csv(rdir / 'top15_per_sample.csv')
        top['sample_id'] = top['sample_id'].astype(str)
        for r in stats.itertuples():
            if r.stage not in ('primary', 'metastatic'):
                continue
            rows = top[top['sample_id'] == str(r.sample_id)].sort_values('rank')
            if rows.empty:
                continue
            fname = f"{r.sample_id}_top15_targets.csv"
            targets = [('all_folder', 'All'), ('all_stage_folder', 'Primary' if r.stage == 'primary' else 'Metastatic')]
            if r.subtype in ('ER+', 'HER2+', 'TNBC'):
                targets.append(('all_subtype_folder', r.subtype))
            for folder, g in targets:
                rows[cols].to_csv(base / folder / g / fname, index=False)
                counts[g] += 1
    print("composite inputs:", counts)
    return base


def patch_gfa(text, folder=''):
    # (1) the units are samples, not patients (124 tumor samples from 109 patients)
    # (2) in the stage folder the two groups are stages, not subtypes
    # (3) presence/absence heatmap: no sample IDs under the columns (100+ overlapping labels were
    #     unreadable); the axis says "Tumor samples (n = N)" instead
    old_ticks = ("    # X-axis labels: first 8 characters of the RENAMED IDs\n"
                 "    short_labels = [pid[:8] for pid in sorted_renamed_ids]\n"
                 "    ax.set_xticks(np.arange(len(short_labels)) + 0.5)\n"
                 "    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)\n")
    new_ticks = ("    # No sample IDs under the columns: 100+ overlapping labels are unreadable\n"
                 "    ax.set_xticks([])\n")
    swaps = [(old_ticks, new_ticks),
             ("ax.set_xlabel('Patient', fontsize=12, fontweight='bold')",
              "ax.set_xlabel(f'Tumor samples (n = {len(sorted_renamed_ids)})', fontsize=12, fontweight='bold')"),
             ("(n={len(sorted_renamed_ids)} patients)", "(n={len(sorted_renamed_ids)} samples)")]
    if folder == 'all_stage_folder':
        swaps.append(('fig.suptitle(f"Per-Subtype Top {top_n}', 'fig.suptitle(f"Per-Stage Top {top_n}'))
    for old, new in swaps:
        assert old in text, f'expected text not found in ORIG_GFA script: {old[:60]!r}'
        text = text.replace(old, new)
    # (4) red labels: known ADC targets instead of the hand-picked fibroblast-marker list
    text, n = re.subn(r"# Known stromal markers to note in plots\nknown_markers = \[[^\]]*\]",
                      "# Known ADC targets (red labels); set by make_original_figures.py\n"
                      f"known_markers = {sorted(KNOWN_ADC_TARGETS)!r}", text)
    assert n == 1, 'known_markers list not found in ORIG_GFA script'
    return text


def run_composite():
    base = build_composite_folders()
    jobs = [('all_folder', ROOT / 'ORIG_GFA' / 'all_folder' / 'gene_frequency_analysis.py'),
            ('all_stage_folder', ROOT / 'ORIG_GFA' / 'all_stage_folder' / 'gene_frequency_analysis_updated.py'),
            ('all_subtype_folder', ROOT / 'ORIG_GFA' / 'all_subtype_folder' / 'gene_frequency_analysis_updated.py')]
    for folder, src in jobs:
        work = base / folder
        script = work / src.name
        script.write_text(patch_gfa(src.read_text(), folder))
        print(f"\n=== composite: running {folder}/{src.name}")
        r = subprocess.run([PY, script.name, '.'], cwd=work, capture_output=True, text=True)
        print('\n'.join(r.stdout.strip().splitlines()[-12:]))
        if r.returncode != 0:
            print(r.stderr[-3000:])
            raise SystemExit(f"{folder} failed")
        out = work / 'gene_frequency_analysis'
        files = sorted(str(p.relative_to(out)) for p in out.rglob('*.png'))
        print(f"  -> {out}: {len(files)} figures")
    eng = ROOT / 'composite' / 'figures_engine'
    f = ROOT / 'composite' / 'figures'
    if f.exists() and not eng.exists():
        shutil.move(str(f), str(eng))


if __name__ == '__main__':
    names = sys.argv[1:] or list(config.DATASETS) + ['composite']
    for n in names:
        if n == 'composite':
            run_composite()
        else:
            run_dataset(n)
