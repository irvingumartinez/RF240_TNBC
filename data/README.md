# RF240 TNBC data directory

This directory contains lightweight derived checkpoints used by the canonical public analysis scripts.

## Included derived checkpoints

- `longitudinal/` — patient-level longitudinal RF240 scores and derived tables supporting Figures 3–4 and Supplementary Figures 4–5.
- `longitudinal/projected_states/` — projected-state composition and RF240 enrichment tables.
- `longitudinal/figure8_crosscohort/` — final Zhang and Shiao support tables used for Figure 8 and Supplementary Figure 8.
- `tcr/figure5/` — clonotype expansion and RF240 module summaries used for Figure 5.
- `tcr/figure6/` — panel-level persistent-clonotype remodeling tables used for Figure 6.
- `tcr/supplementary_figure6/` — repertoire diversity, public/private clonotype, and expansion-category support tables.
- `bassez/figure7/` — audited Bassez Figure 7 plotting and support tables.
- `bassez/supplementary_figure7/` — audited Bassez patient-level robustness and transition-context tables.

## Baseline objects retained locally but excluded from GitHub

The following files may exist locally but are intentionally excluded from normal Git tracking:

- `baseline/GSE246613_Tcells_Base_only.h5ad`
- `baseline/sce_tbase_curated_reanalysis.rds`
- `baseline/meta_reanalysis_baseline_Tcells.csv`
- `baseline/metadata_Tcells_Base_RF240.csv`

The `.h5ad` is a large source object and the remaining files are large curated or derived checkpoints. Their provenance and analytical role are described in `../documentation/DATA_PROVENANCE.md`.

These files should not be committed to GitHub unless a deliberate redistribution decision is made. Large derived checkpoints may instead be deposited in an archival repository such as Zenodo at publication time.

## Reproduction strategy

Scripts for Figures 1–2 and Supplementary Figures 2–3 depend on the curated baseline object or associated baseline checkpoints. Scripts for later figures primarily use the smaller derived tables included here.

No source dataset should be redistributed solely because it is present locally; original licensing, repository terms, and publication policy should be respected.
