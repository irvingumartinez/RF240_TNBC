# RF240 TNBC code compendium

This repository contains the analysis code, metadata, and lightweight derived checkpoints supporting the manuscript:

**Biological deconstruction of an RF-selected T-cell signature reveals state- and clonotype-dependent organization during immunotherapy in triple-negative breast cancer**

## Scope

The study deconstructs an RF-selected 240-gene T-cell signature (RF240) associated with immunotherapy response in triple-negative breast cancer and examines its organization across transcriptional states, longitudinal treatment, TCR clonotypes, and persistent-clonotype remodeling.

The repository is organized around figure-level reproducibility. Final public scripts reproduce manuscript panels from audited derived checkpoints, while upstream provenance and processing relationships are documented separately.

## Repository structure

- `scripts/` — canonical public analysis scripts (`00_setup.R` through `15_SuppFigure8_crosscohort_sensitivity.R`).
- `metadata/` — RF240 gene list and curated baseline cluster identity map.
- `data/` — lightweight derived checkpoints used by the public figure scripts.
- `results/` — output directories for regenerated figures, supplementary figures, and support tables.
- `documentation/` — analysis workflow, data provenance, figure-to-script map, and pre-GitHub file classification.

## Reproducibility model

This compendium uses two complementary reproducibility layers:

1. **Figure-level reproducibility** — each canonical script reproduces a final manuscript figure or supplementary figure from an audited derived checkpoint.
2. **Upstream traceability** — the origin of those checkpoints, source datasets, processing steps, and historical analysis scripts is documented in `documentation/DATA_PROVENANCE.md` and `documentation/ANALYSIS_WORKFLOW.md`.

Heavy source objects and large baseline checkpoints are intentionally not tracked in GitHub. Their filenames and roles are documented in `data/README.md` and `documentation/DATA_PROVENANCE.md`.

## Execution

Start by sourcing:

```r
source("scripts/00_setup.R")
```

Then run the desired figure script from `scripts/`. The complete figure-to-script mapping is provided in `documentation/FIGURE_SCRIPT_MAP.md`.

## Statistical interpretation

Patient-level inference is treated as primary where repeated clonotypes are nested within individuals. Clone-level analyses are retained as exploratory or contextual where appropriate. Bassez expansion-associated (`E`) and non-expansion-associated (`NE`) groups are not equivalent to clinical responders and nonresponders.

## Data availability

Large raw or curated baseline objects are excluded from normal Git tracking. Public-source datasets should be obtained from their original repositories where applicable; derived checkpoints required for figure-level reproduction are included when sufficiently small and appropriate for redistribution.

See `data/README.md` and `documentation/DATA_PROVENANCE.md` for details.
