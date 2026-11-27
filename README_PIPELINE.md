# Pipeline details

All six scRNA-seq datasets are analyzed by one shared engine (`common/`), driven
by one configuration entry per dataset (`common/config.py`).

```
python run_pipeline.py all              # steps 2–6 for every dataset (~10 min total)
python run_pipeline.py GSE161529        # one dataset
python run_pipeline.py GSE305078 --steps 4,5,6
python composite_scoring.py             # cross-dataset composite score
python composite_scoring.py --legacy ORIG_GFA   # recompute the original hand-made analysis
```

Each run archives the previous `results/` contents into `results/_archive_<date>/`
and writes `results/pipeline_run.log`.

## What a run produces (identical names in every dataset)

| file | content |
|---|---|
| `celltype_mapping.json` | stromal/tumor definition actually used, cell counts |
| `surface_genes.txt` | the GO-derived surface universe for this dataset |
| `global_stromal_targets.csv` | all surface genes with ≥10 % stromal detection, ranked |
| `<subtype>_stromal_targets.csv`, `subtype_comparison.csv` | per-subtype scoring (≥2 samples & ≥50 stromal cells) |
| `per_sample_stromal_targets.csv`, `sample_statistics.csv` | per-sample scores of the top-200 global genes; every sample listed with `analyzed` flag |
| `individual_samples/<sample>_top15_targets.csv` | per-sample top 15 |
| `primary_stromal_targets.csv`, `metastatic_stromal_targets.csv`, `primary_vs_metastatic_comparison.csv` | stage aggregation (sample level) |
| `top15_per_sample.csv`, `top15_gene_frequencies.csv`, `top15_summary.txt` | inputs to the composite score |
| `figures/fig1…fig5*.pdf`, `figures/table1_top_targets.csv` | figures |

`composite/` holds the pooled analysis: one folder per group (All, Primary,
Metastatic, ERpos, HER2pos, TNBC) with `gene_frequency_summary.csv` and
`top25_composite.txt` (the list to feed the TCGA/METABRIC/MET500 scripts),
plus `sample_counts.csv` (the *n* to quote) and `samples_used.csv`.

## Definitions (identical across datasets)

* **stromal** = fibroblast/CAF lineage only (pericytes included where the source
  annotation lumps them in). Endothelial cells are never stromal.
* **tumor** = malignant / all epithelial clusters.
* **surface genes** = genes annotated in Gene Ontology to the plasma membrane and
  to the cell surface, minus secreted/extracellular-matrix terms. Whitelisted
  surface families (CDH\*, ITGA\*, ITGB\*, CD##, PDGFRA/B, FAP, THY1, PDPN, …) are
  kept whenever they are membrane-anchored, even if a shed or soluble form is also
  annotated as extracellular (see `common/surface.py`).
* **expression** = log1p(CP10k). Datasets whose h5ad stores raw UMI counts
  (GSE167036, GSE176078, GSE306201) are normalized on load.
* **score** = `log2FC × stromal_detection × tumor_penalty × detection_specificity`
  on log1p(CP10k) data, with `detection_specificity = min(stromal_det / all_det, 10)`
  and set to 0 when the gene is detected in < 10 stromal cells (see `common/scoring.py`).
* per-sample analysis: samples with ≥ 20 stromal cells; top-200 global candidates.
* **stage** is assigned per sample, so a patient's primary tumor and lymph-node
  samples are analyzed separately.
* composite score = `frequency_pct / mean_rank` over per-sample top-15 lists
  (same formula as the original `ORIG_GFA` script; tumor samples only).
