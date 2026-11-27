"""
Per-dataset configuration.

Every dataset is described by the same small set of fields so that the shared
engine (loading/scoring/steps) treats all six identically.  Anything that is
dataset-specific (which obs column holds the cell type, which labels are
stromal, how stage/subtype are derived) lives here and nowhere else.

Definitions used throughout
---------------------------
stromal  : fibroblast / CAF lineage ONLY (pericytes included where the source
           annotation lumps them with fibroblasts).  Endothelial cells are NOT
           stromal in any dataset (previously GSE167036 and GSE306201 included
           them, which pushed endothelial genes such as CDH5/CLDN5/PLVAP into
           the "metastatic" target list).
tumor   : malignant / cancer epithelial cells (all epithelial clusters in
           datasets without a malignant-vs-normal epithelial call).
stage    : primary | metastatic | normal | preneoplastic
subtype  : ER+ | HER2+ | TNBC | PR+ | Normal | BRCA1-preneoplastic
"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Scoring constants (single source of truth; see scoring.py for the formula)
# ---------------------------------------------------------------------------
PSEUDOCOUNT = 0.001
MIN_STROMAL_DETECTED_CELLS = 10     # gene must be detected in >=10 stromal cells
SPECIFICITY_CAP = 10.0              # cap on stromal_detection / all_detection
MIN_TUMOR_CELLS_FOR_PENALTY = 10    # tumor penalty only evaluated with >=10 tumor cells
MIN_STROMAL_DETECTION = 0.10        # reporting filter: >=10 % of stromal cells
MIN_STROMAL_CELLS_PER_SAMPLE = 20   # per-sample analysis requires >=20 stromal cells
N_GLOBAL_CANDIDATES = 200           # per-sample analysis scores the top-200 global genes
TOP_N_PER_SAMPLE = 15
MIN_STROMAL_CELLS_PER_SUBTYPE = 50
MIN_SAMPLES_PER_SUBTYPE = 2

# Known surface-expressed CAF markers used only for sanity checks.
# (COL1A1 / POSTN are secreted ECM proteins and are excluded from the surface
#  universe by design, so they must NOT be in this list.)
KNOWN_STROMAL_SURFACE_MARKERS = ['FAP', 'PDPN', 'THY1', 'PDGFRA', 'PDGFRB',
                                 'LRRC15', 'CDH11', 'ITGA11', 'MXRA8', 'ANTXR1']
EPITHELIAL_MARKERS = ['EPCAM', 'KRT19', 'KRT18', 'CDH1', 'MUC1', 'NECTIN4', 'ERBB2']


def _gse180286_subtype(obs):
    # From GSE180286 series matrix "molecular typing": P1 TNBC, P2 luminal B,
    # P3 Her-2, P4 TNBC, P5 Her-2.  Patient ids in the h5ad are A2019..E2020.
    m = {'A2019': 'TNBC', 'B2019': 'ER+', 'C2020': 'HER2+',
         'D2020': 'TNBC', 'E2020': 'HER2+'}
    return obs['patient_id'].astype(str).map(m).fillna('Unknown')


def _gse167036_subtype(obs):
    er = obs['ER_statu'].astype(str).str.lower().str.startswith('pos')
    her2 = obs['Her2_statu'].astype(str).str.lower().str.startswith('pos')
    out = er.map({True: 'ER+', False: 'Unknown'})
    out[her2] = 'HER2+'
    return out


def _gse306201_subtype(obs):
    # All patients in this cohort are ER+ (MajorSubtype == 'ER').
    return obs['MajorSubtype'].astype(str).map({'ER': 'ER+'}).fillna('Unknown')


DATASETS = {
    # ------------------------------------------------------------------ #
    'GSE161529': dict(
        h5ad='GSE161529/data/dataset.h5ad',
        results='GSE161529/results',
        celltype_col='celltype_annotation',
        stromal_labels=['Fibroblast'],
        tumor_labels=['Epithelial'],
        sample_col='sample_id',
        patient_col='patient_id',
        # obs columns written by fix_metadata_gse161529.py (from the GEO series matrix)
        stage_col='stage',
        subtype_col='cancer_subtype',
        description='Pal et al. 2021 – 69 samples: normal, BRCA1 pre-neoplastic, '
                    'ER+/PR+/HER2+/TNBC primaries and involved lymph nodes',
    ),
    # ------------------------------------------------------------------ #
    'GSE167036': dict(
        h5ad='GSE167036/data/dataset.h5ad',
        results='GSE167036/results',
        celltype_col='major_celltype_name',
        stromal_labels=['CAFs'],
        tumor_labels=['Epithelial Cells'],
        sample_col='orig.ident',            # CA1..CA8 (tumor), LN1..LN8 (lymph node)
        patient_col='patient_id',
        stage_map=('sample_type', {'Tumor': 'primary', 'Lymph Node': 'metastatic'}),
        subtype_fn=_gse167036_subtype,
        description='8 patients, matched primary tumor + lymph-node metastasis',
    ),
    # ------------------------------------------------------------------ #
    'GSE176078': dict(
        h5ad='GSE176078/data/dataset.h5ad',
        results='GSE176078/results',
        celltype_col='celltype_major',
        stromal_labels=['CAFs'],
        tumor_labels=['Cancer Epithelial'],
        sample_col='orig.ident',
        patient_col='orig.ident',
        stage_const='primary',
        subtype_col='subtype',
        description='Wu et al. 2021 – 26 primary tumors (ER+/HER2+/TNBC)',
    ),
    # ------------------------------------------------------------------ #
    'GSE180286': dict(
        h5ad='GSE180286/data/dataset.h5ad',
        results='GSE180286/results',
        # preprocess.py had no 'Endothelial' category, so the endothelial
        # cluster (leiden 17: PECAM1/VWF/CLDN5/PLVAP high, COL1A1 low) was
        # labeled 'Fibroblast'.  Use the leiden clusters directly:
        #   fibroblast/CAF: 5, 20, 22 ; pericyte-like (ACTA2/RGS5): 21
        #   epithelial    : 4, 9, 10, 13, 15, 18 (== celltype_annotation 'Epithelial')
        celltype_col='leiden',
        stromal_labels=['5', '20', '21', '22'],
        tumor_labels=['4', '9', '10', '13', '15', '18'],
        sample_col='sample_id',
        patient_col='patient_id',
        stage_map=('sample_type', {'primary': 'primary', 'metastatic': 'metastatic'}),
        subtype_fn=_gse180286_subtype,
        description='5 patients, primary tumor + 2 lymph nodes each',
    ),
    # ------------------------------------------------------------------ #
    'GSE305078': dict(
        h5ad='GSE305078/data/GSE305078_dataset.h5ad',
        results='GSE305078/results',
        celltype_col='leiden',
        # From results/celltype_scores.csv (marker-based): clusters with
        # >=50 % stromal-marker share and stromal score >= 0.5.
        # Cluster 18 (87.6 % stromal) was previously dropped only because it
        # fell 125 cells short of a 1 %-of-dataset size threshold.
        stromal_labels=['6', '15', '18', '25'],
        # All epithelial clusters (>=50 % epithelial-marker share); the old
        # pipeline used only cluster 8, the LEAST epithelial of these.
        tumor_labels=['2', '5', '8', '9', '11', '12', '14', '17', '19', '20'],
        sample_col='sample_id',
        patient_col='patient_id',
        # GEO: "early-stage TNBC" vs "advanced TNBC" pre-treatment tumor biopsies.
        # 'advanced' is mapped to the metastatic arm for continuity with the
        # previous analysis; see stage_detail in sample_statistics.csv.
        stage_map=('sample_type', {'early_stage': 'primary', 'advanced': 'metastatic'}),
        stage_detail_map=('sample_type', {'early_stage': 'early-stage TNBC',
                                          'advanced': 'advanced TNBC (pre-treatment biopsy)'}),
        subtype_const='TNBC',
        description='32 TNBC biopsies (26 early-stage, 6 advanced)',
    ),
    # ------------------------------------------------------------------ #
    'GSE306201': dict(
        h5ad='GSE306201/data/dataset.h5ad',
        results='GSE306201/results',
        celltype_col='major_celltype',
        stromal_labels=['Fibroblast'],
        tumor_labels=['Malignant'],
        sample_col='Sample',
        patient_col='Sample',
        stage_map=('Disease', {'Primary': 'primary', 'Metastasis': 'metastatic'}),
        stage_detail_col='MetastaticTarget',
        subtype_fn=_gse306201_subtype,
        description='ER+ cohort – 12 primaries, 11 distant metastases',
    ),
}


def get(name):
    if name not in DATASETS:
        raise KeyError(f"Unknown dataset '{name}'. Known: {', '.join(DATASETS)}")
    cfg = dict(DATASETS[name])
    cfg['name'] = name
    cfg['h5ad'] = PROJECT_ROOT / cfg['h5ad']
    cfg['results'] = PROJECT_ROOT / cfg['results']
    return cfg
