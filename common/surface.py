"""
Cell-surface gene universe from GO cellular-component annotations (GOA human).

Set algebra (unchanged from the original scripts except for the rescue rule):
  A = membrane-anchored   : GO:0005886 (plasma membrane) or GO:0005887 (integral component of PM)
  B = accessible surface  : GO:0009897 (external side of PM) or GO:0009986 (cell surface)
  E = secreted / ECM      : GO:0005576, GO:0005615, GO:0031012
  core    = (A & B) - E
  rescue  = whitelisted surface families (CDH*, ITGA*, ITGB*, CD##, PDGFRA/B, FAP, THY1, PDPN ...)
  surface = core | (rescue & A)

Rescue rule change: the original code used `rescue & (A - E)`, i.e. the
exclusion list was applied BEFORE the whitelist, so FAP, THY1, CDH11, LRRC15, CD34
etc. (membrane proteins that also carry an "extracellular region" annotation
because they shed/soluble forms exist) were silently dropped in every dataset.
The whitelist now only requires membrane anchoring (A).

Cache: one GAF copy + one parsed-term cache shared by all datasets under
common/cache/.  The per-dataset surface set is keyed on a hash of the dataset's
gene list and RULE_VERSION so a rule change can never reuse a stale set.
"""
import gzip
import hashlib
import json
import pickle
import re
import shutil
from datetime import datetime
from pathlib import Path

from . import config

RULE_VERSION = 'v2-rescue-requires-A-only'
CACHE_DIR = config.PROJECT_ROOT / 'common' / 'cache'
GAF_FILE = CACHE_DIR / 'goa_human.gaf.gz'
GAF_URL = 'https://current.geneontology.org/annotations/goa_human.gaf.gz'

MEMBRANE_ANCHORED = {'GO:0005886', 'GO:0005887'}
ACCESSIBLE_SURFACE = {'GO:0009897', 'GO:0009986'}
EXCLUDE_SECRETED_ECM = {'GO:0005576', 'GO:0005615', 'GO:0031012'}

_PREFIXES = ('CDH', 'ITGA', 'ITGB', 'TSPAN', 'SDC', 'GPC', 'CLDN')
_REGEXES = [re.compile(p) for p in (
    r'^(ICAM|VCAM)\d*$',
    r'^(EGFR|ERBB\d)$',
    r'^(FGFR\d)$',
    r'^(PDGFRA|PDGFRB)$',
    r'^(KDR|FLT1|FLT4)$',
    r'^(EPHA\d+|EPHB\d+)$',
    r'^(CXCR\d+|CCR\d+|ACKR\d+)$',
    r'^CD\d+[A-Z0-9]*$',
)]
_SINGLE = {'EPCAM', 'PROM1', 'PECAM1', 'ENG', 'MCAM', 'NT5E', 'PDPN', 'THY1', 'FAP',
           'CD44', 'CD47', 'CD55', 'CD59', 'LRRC15'}


def _rescue_family(genes):
    return {g for g in genes
            if g in _SINGLE or g.startswith(_PREFIXES) or any(rx.match(g) for rx in _REGEXES)}


def _ensure_gaf(log):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if GAF_FILE.exists():
        return
    # reuse a copy already downloaded by the old per-dataset scripts
    for cand in sorted(config.PROJECT_ROOT.glob('GSE*/results/goa_human.gaf.gz')):
        log(f"  copying existing GAF from {cand}")
        shutil.copy(cand, GAF_FILE)
        return
    import requests
    log(f"  downloading {GAF_URL}")
    r = requests.get(GAF_URL, stream=True, timeout=300)
    r.raise_for_status()
    tmp = GAF_FILE.with_suffix('.tmp')
    with open(tmp, 'wb') as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    tmp.replace(GAF_FILE)


def _gene2cc(log):
    """gene symbol -> set(GO CC ids), parsed once and cached."""
    cache = CACHE_DIR / 'gene2cc.pkl'
    if cache.exists():
        with open(cache, 'rb') as f:
            return pickle.load(f)
    _ensure_gaf(log)
    log('  parsing GOA cellular-component annotations')
    g2c = {}
    with gzip.open(GAF_FILE, 'rt') as f:
        for line in f:
            if line.startswith('!'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 15 or cols[8] != 'C' or 'NOT' in (cols[3] or ''):
                continue
            g2c.setdefault(cols[2], set()).add(cols[4])
    with open(cache, 'wb') as f:
        pickle.dump(g2c, f, protocol=pickle.HIGHEST_PROTOCOL)
    return g2c


def get_surface_genes(var_names, results_dir, log=print, force=False):
    genes = sorted(set(map(str, var_names)))
    key = hashlib.sha1(('\n'.join(genes) + RULE_VERSION).encode()).hexdigest()[:12]
    results_dir = Path(results_dir)
    cache = results_dir / f'surface_genes_{key}.json'
    if cache.exists() and not force:
        with open(cache) as f:
            d = json.load(f)
        log(f"  surface universe from cache ({len(d['surface_genes'])} genes, rule {RULE_VERSION})")
        return d['surface_genes'], d['meta']

    g2c = _gene2cc(log)
    gs = set(genes)
    in_goa = gs & set(g2c)
    A = {g for g in in_goa if g2c[g] & MEMBRANE_ANCHORED}
    B = {g for g in in_goa if g2c[g] & ACCESSIBLE_SURFACE}
    E = {g for g in in_goa if g2c[g] & EXCLUDE_SECRETED_ECM}
    core = (A & B) - E
    rescue = _rescue_family(gs) & A
    surface = sorted(core | rescue)
    meta = {
        'rule_version': RULE_VERSION,
        'n_dataset_genes': len(genes),
        'n_with_goa_cc': len(in_goa),
        'n_membrane_anchored_A': len(A),
        'n_accessible_B': len(B),
        'n_excluded_E': len(E),
        'n_core': len(core),
        'n_rescued_total': len(rescue),
        'n_rescued_beyond_core': len(rescue - core),
        'n_surface': len(surface),
        'rescued_beyond_core': sorted(rescue - core),
        'timestamp': datetime.now().isoformat(timespec='seconds'),
    }
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(cache, 'w') as f:
        json.dump({'surface_genes': surface, 'meta': meta}, f, indent=1)
    with open(results_dir / 'surface_genes.txt', 'w') as f:
        f.write('\n'.join(surface) + '\n')
    log(f"  surface universe built: {len(surface)} genes "
        f"(core {len(core)} + rescued {len(rescue - core)}); rule {RULE_VERSION}")
    missing = [g for g in config.KNOWN_STROMAL_SURFACE_MARKERS if g in gs and g not in surface]
    if missing:
        log(f"  WARNING known surface CAF markers not in universe: {missing}")
    return surface, meta
