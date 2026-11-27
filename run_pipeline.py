#!/usr/bin/env python
"""
Run the shared stromal-target pipeline for one or more datasets.

    python run_pipeline.py GSE180286                 # steps 2-6
    python run_pipeline.py GSE161529 --steps 2,3     # subset
    python run_pipeline.py all                       # every dataset
    python run_pipeline.py all --no-archive          # keep existing results in place

Each run writes results/<dataset>/pipeline_run.log and archives any previous
results into results/_archive_<date>/ first (unless --no-archive).
Steps 2 and 3 need the h5ad in memory; 4-6 work from the CSVs.
"""
import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import config, loading, steps, figures  # noqa: E402


class Logger:
    def __init__(self, path):
        self.f = open(path, 'a')

    def __call__(self, msg=''):
        print(msg, flush=True)
        self.f.write(msg + '\n')
        self.f.flush()


def run(name, step_list, archive=True, force_surface=False):
    cfg = config.get(name)
    cfg['results'].mkdir(parents=True, exist_ok=True)
    if archive and min(step_list) == 2:
        steps.archive_old_results(cfg, print)
    log = Logger(cfg['results'] / 'pipeline_run.log')
    log(f"\n{'#' * 80}\n# {name}  {datetime.now():%Y-%m-%d %H:%M}  steps {step_list}\n# {cfg['description']}\n{'#' * 80}")
    t0 = time.time()
    adata = None
    if any(s in step_list for s in (2, 3)):
        adata = loading.load_dataset(cfg, log)
        loading.write_celltype_mapping(cfg, adata)
    if 2 in step_list:
        steps.step2_global(cfg, adata, log, force_surface=force_surface)
    if 3 in step_list:
        steps.step3_samples(cfg, adata, log)
    del adata
    if 4 in step_list:
        steps.step4_stage(cfg, log)
    if 5 in step_list:
        figures.step5_figures(cfg, log)
    if 6 in step_list:
        steps.step6_top15(cfg, log)
    log(f"\n{name} done in {(time.time() - t0) / 60:.1f} min")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('datasets', nargs='+')
    ap.add_argument('--steps', default='2,3,4,5,6')
    ap.add_argument('--no-archive', action='store_true')
    ap.add_argument('--force-surface', action='store_true', help='rebuild the GO surface set')
    a = ap.parse_args()
    names = list(config.DATASETS) if a.datasets == ['all'] else a.datasets
    step_list = sorted(int(s) for s in a.steps.split(','))
    for n in names:
        run(n, step_list, archive=not a.no_archive, force_surface=a.force_surface)
