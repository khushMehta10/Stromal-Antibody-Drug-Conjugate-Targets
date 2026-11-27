#!/usr/bin/env python
"""Run the shared pipeline for GSE180286 (steps 2-6 by default).

    python run_GSE180286.py              # all steps
    python run_GSE180286.py --steps 4,5,6

Dataset definitions live in common/config.py; see README_PIPELINE.md.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_pipeline import run  # noqa: E402

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', default='2,3,4,5,6')
    ap.add_argument('--no-archive', action='store_true')
    a = ap.parse_args()
    run('GSE180286', sorted(int(s) for s in a.steps.split(',')), archive=not a.no_archive)
