# Analysis workflow

## Reproducibility strategy

The RF240 TNBC compendium is organized around figure-level reproducibility from audited checkpoints rather than mandatory re-execution of every upstream heavy analysis step. Upstream provenance is documented in `DATA_PROVENANCE.md`.

## Recommended execution order

1. `00_setup.R` — initialize repository-relative paths and shared project settings.
2. `01_Figure1_baseline_Tcell_architecture.R` — baseline T-cell architecture and RF240 distribution.
3. `02_SuppFigure2_RF240_spatial_support.R` — extended baseline RF240 visualization and ordered gene support.
4. `03_SuppFigure3_functional_enrichment.R` — state/module functional enrichment from the curated baseline object.
5. `04_Figure2_RF240_coexpression_architecture.R` — patient-level RF240 coexpression architecture.
6. `05_Figure3_longitudinal_RF240_dynamics.R` — longitudinal global and module-level RF240 dynamics.
7. `06_SuppFigure4_patient_trajectories.R` — patient-level global and module RF240 trajectories.
8. `07_Figure4_state_resolved_RF240.R` — projected T-cell states and RF240 enrichment across treatment.
9. `08_SuppFigure5_state_resolved_support.R` — extended state-resolved longitudinal context.
10. `09_Figure5_clonotype_RF240_organization.R` — state-resolved clonal expansion and RF240 organization.
11. `10_SuppFigure6_repertoire_structure.R` — public/private clonotype structure, diversity, and RF240 module organization.
12. `11_Figure6_persistent_clonotype_remodeling.R` — patient- and clonotype-level RF240 remodeling across longitudinal clonotype fates.
13. `12_Figure7_Bassez_contextualization.R` — phenotype-state and RF240 remodeling in persistent clonotypes from the Bassez cohort.
14. `13_SuppFigure7_Bassez_robustness.R` — patient-level Bassez robustness, expansion context, and transition summaries.
15. `14_Figure8_crosscohort_RF240_remodeling.R` — cross-cohort RF240 remodeling in persistent clonotypes.
16. `15_SuppFigure8_crosscohort_sensitivity.R` — Zhang/Shiao sensitivity analyses and negative-control context.

## Dependency logic

The scripts are numbered in manuscript order, but not every script requires all previous scripts to be rerun. Many later scripts read audited derived tables committed under `data/`.

The main dependency that should be preserved is the baseline ordering support used across the early figures: Supplementary Figure 2 provides the corrected RF240 gene ordering used downstream in baseline-related visualization.

## Analytical hierarchy

- Figures 1–2 establish baseline T-cell state architecture and RF240 coexpression organization.
- Figures 3–4 place RF240 in longitudinal and projected-state treatment context.
- Figure 5 links RF240 organization to TCR clonotype structure and expansion.
- Figure 6 examines RF240 remodeling across persistent clonotype fates, with patient-level inference treated as primary.
- Figure 7 contextualizes persistent-clonotype remodeling in the Bassez cohort.
- Figure 8 summarizes Zhang and Shiao cross-cohort contextualization and robustness.

## Statistical guardrails

- Patient-level inference is primary when clonotypes are nested within patients.
- Clone-level analyses are exploratory/contextual unless explicitly stated otherwise.
- Bassez `E` and `NE` denote expansion-associated and non-expansion-associated clonotypes, not clinical responders and nonresponders.
- Cross-cohort analyses provide contextual support and should not be described as independent validation when source relationships or sample structure do not justify that claim.
- Correlation and network analyses are descriptive and are not interpreted causally.

## Output organization

Each script writes regenerated plots and support tables under `results/figures/`, `results/supplementary_figures/`, or `results/tables/` using repository-relative paths.

For exact figure-to-script correspondence, see `FIGURE_SCRIPT_MAP.md`.
