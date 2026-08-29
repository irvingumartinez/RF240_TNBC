# Data provenance

This document records the analytical provenance of the RF240 TNBC manuscript.

The repository is designed around audited derived checkpoints. It does not require every reviewer or user to reproduce all single-cell preprocessing from raw files before regenerating figures. Upstream provenance remains documented, while figure-generation scripts use the smallest audited checkpoint that preserves the reported analysis.

## Reproducibility model

The repository has two reproducibility layers.

### Layer 1: figure-level reproducibility
Figure scripts regenerate the reported statistics and visual outputs from audited derived inputs.

### Layer 2: upstream analytical traceability
Upstream scripts and documentation describe how those derived inputs were produced from the original public datasets, including filtering, normalization, scoring, state assignment, clonotype definition, and longitudinal summaries.

This distinction avoids unnecessary reprocessing of large single-cell objects while preserving a transparent analytical chain.

## Dataset and checkpoint provenance

### Baseline TNBC T-cell reanalysis

Primary baseline single-cell input:
- `GSE246613_Tcells_Base_only.h5ad`

Audited preprocessing pathway:
- cell QC: 300–5000 detected genes;
- total counts <= 20,000;
- mitochondrial percentage <= 10%;
- APC-like residual score based on `CD74`, `HLA-DRA`, `HLA-DPA1`, `IFI30`, `APOE`, `CST3`, with removal above the 99th percentile;
- size-factor normalization and log-normalization;
- HVG selection: 3,000 genes;
- exact PCA: 30 PCs;
- k-means clustering: 8 centers using the first 20 PCs, `nstart = 20`, seed 123;
- UMAP: `n_neighbors = 20`, `min_dist = 0.3`.

Audited checkpoint:
- `sce_tbase_curated_reanalysis.rds`
- `meta_reanalysis_baseline_Tcells.csv`

Public fixed cluster map:
- C1: IFN-high cytotoxic
- C2: Activated CD4/Treg-like
- C3: APC-like residual
- C4: Ambiguous residual
- C5: Cytotoxic CD8 effector-like
- C6: Activated helper/memory IEG-like
- C7: gdT/NK-like terminal cytotoxic
- C8: Helper/memory basal-like

RF240 gene list:
- public filename: `metadata/RF240_genes.csv`

Important audited corrections:
- Figure 1E uses the final dominant-cluster ordering correction.
- Figure 1F uses the true overall cluster proportion, `n_total_cells / sum(n_total_cells)`, not the historical sum of within-response proportions.

### Baseline RF240 coexpression and enrichment

Inputs are derived from the curated baseline T-cell checkpoint and RF240 gene list.

Outputs support:
- Figure 2 coexpression architecture;
- Supplementary Figure 2 extended RF240 spatial/state support;
- Supplementary Figure 3 functional enrichment.

All public versions should use the audited curated checkpoint and final RF240 ordering rules.

### Longitudinal RF240 scoring

Audited longitudinal score products include:
- patient-time global RF240 score;
- patient-time module scores;
- patient-time gene-level matrix where required by upstream logic;
- projected state assignments;
- state-resolved RF240 summaries.

Public Supplementary Figure 4 uses only patient-level global/module trajectories and should not duplicate Figure 3 tests.

### Longitudinal projected-state analysis

Figure 4 and Supplementary Figure 5 use the audited projected-state mapping and RF240 state summaries.

The public state-identity map should be stored under:
- `metadata/projected_state_identity_map.csv`

The public script should preserve the final state mapping and final Figure 4 visual range/label decisions.

### TCR clonotype and repertoire analyses

Figure 5 / Supplementary Figure 6 use audited TCR-linked cell and clonotype summaries.

Important Supplementary Figure 6 correction:
- public-clonotype fractions must use the exact intended denominator columns:
  - `prop_clonotypes`
  - `prop_tcr_cells`
  - `prop_cells_cluster`
- absent public rows are represented as zero, not missing.

Diversity metrics should be selected explicitly by name rather than by positional selection.

### Persistent clonotype remodeling

Figure 6 uses audited persistent-clonotype fate and RF240 remodeling summaries.

Primary inputs include:
- `clonotype_fate_type_summary.csv`
- `RF240_clone_level_module_score_deltas_by_fate.csv`
- `RF240_transition_summary_by_module.csv`

Primary inference is patient-level.
Clone-level results are exploratory/contextual.

### Bassez cohort contextualization

The Bassez analyses are contextual/supportive and should not be described as independent validation.

Final audited clone source:
- `56F_C_Bassez_clone_delta_state_change_finalQC_min100cells_min10persistent.csv.gz`

Primary QC:
- at least 100 T cells at both timepoints;
- at least 10 persistent clonotypes.

Final QC universe:
- 1,506 persistent clonotypes;
- 21 patients.

Expansion terminology:
- E = expansion-associated;
- NE = non-expansion-associated.

Key audited patient-level result:
- E: n = 11, median DeltaRF240 = 0.0263
- NE: n = 10, median DeltaRF240 = -0.0187
- nominal Wilcoxon p = 0.0317
- historical adjusted q = 0.0952
- bootstrap CI for E-minus-NE difference: 0.000448 to 0.0837

Important correction:
- `63D_3b` corrects a descriptive-summary bug in `63D_3`; the primary inference itself was unaffected.
- Public code should contain only the corrected implementation.

Transition summaries in the final public Supplementary Figure 7 should use the final-QC universe consistently.

### Zhang contextual analysis

Patient-level RF240 remodeling:
- NR P017: -0.104
- NR P005: -0.0565
- NR P002: -0.0140
- R P012: +0.0168
- R P019: +0.0939

Interpretation:
- directional separation is descriptive;
- exact Wilcoxon p = 0.20;
- descriptive bootstrap interval remains positive across patient resampling;
- leave-one-out and minimum-clone threshold sensitivity retain the same direction.

This is contextual support, not independent validation.

### Shiao contextual analysis

Final patient-level QC:
- Base to PD1: n = 25 (5 NR, 20 R)
- Base to RTPD1: n = 20 (4 NR, 16 R)
- PD1 to RTPD1: n = 22 (5 NR, 17 R)

No response comparison remains significant after adjustment.
A weak late nominal signal is retained only as contextual evidence.
No consistent state-switch, phenotype-coupling, or expansion effect is claimed.

### Cross-dataset synthesis

Figure 8 combines Zhang and Shiao patient-level longitudinal context.

The figure is designed as cross-dataset contextualization, not as an independent validation analysis.

Public wording should avoid:
- `validation`;
- `replication`;
- claims of independent predictive confirmation.

## Statistical-unit policy

- Patient-level inference is primary whenever patients are the independent biological units.
- Clone-level tests are exploratory/contextual unless explicitly stated otherwise.
- Bassez E/NE labels indicate expansion association, not treatment response.
- Correlation networks are descriptive and should not be interpreted causally.
- `RF240 remodeling` is reserved primarily for clonotype-level longitudinal DeltaRF240 analyses.

## Data-distribution policy

Every input should be classified before public release as one of:

1. redistributable;
2. derived and shareable;
3. must be downloaded from the original public repository;
4. too large for GitHub and therefore hosted externally or reconstructed;
5. restricted or not redistributable.

Large files such as `.h5ad` should not be committed directly to GitHub.

The current `.gitignore` excludes large raw single-cell formats.

## Repository-cleaning rule

Public scripts should preserve the audited analytical logic while removing:
- obsolete branches;
- debugging code;
- internal provenance diagnostics;
- historical version suffixes;
- temporary paths;
- duplicated plotting experiments.

Corrections discovered during audit should be integrated directly into the single canonical public script rather than represented as chains of patch scripts.
