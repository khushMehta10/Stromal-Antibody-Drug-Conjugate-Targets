"""
Shared stromal ADC target-discovery engine used by all six scRNA-seq datasets.

Modules
-------
config   : per-dataset definitions (file paths, stromal/tumor labels, stage & subtype columns)
loading  : h5ad loading, normalization check, standardized obs columns
surface  : GO-based cell-surface gene universe (shared cache)
scoring  : vectorized stromal-exclusivity score (single source of truth for the formula)
steps    : step2 (global/subtype), step3 (per-sample), step4 (stage), step6 (top-15)
figures  : step5 figures
"""
