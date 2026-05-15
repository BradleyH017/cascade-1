# CASCADE

[![R >= 3.5.0](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**C**omprehensive **A**nalysis **S**uite for **C**ell-type specificity and **A**ccessibility-**D**riven QTL **E**ffects

## Overview

CASCADE (Comprehensive Analysis Suite for Cell-type specificity and Accessibility-Driven QTL Effects) is an R package for systematically characterizing cell type specificity and regulatory mechanisms of molecular QTLs (molQTLs) from multiome analyses. It implements a hierarchical framework that operates at two levels:

1. **Gene- or peak-level analysis** classifies eGenes and caPeaks by their cell type specificity patterns, incorporating power-aware assessment through multivariate adaptive shrinkage ([mash](https://github.com/stephenslab/mashr)). Specificity is defined hierarchically over cell lineages, and local false sign rates (LFSR) are used to identify likely shared but underpowered features whose apparent specificity reflects power rather than biology.
2. **Variant-level analysis** classifies individual fine-mapped variants by their regulatory mechanism (chromatin-to-expression cascades) and cell type specificity. Twenty-five mutually exclusive QTL patterns are grouped into eight mechanism categories, ordered by how directly the variant acts on its own regulatory sequence.

## Installation

```r
remotes::install_github("mkanai/cascade")
```

### Dependencies

**Required:** `cli`, `data.table (>= 1.12.0)`, `digest`, `jsonlite`, `memoise`, `parallel`, `qs2`, `R.utils`, `Rcpp`

**Suggested:** `GenomicRanges`, `IRanges`, `S4Vectors` (peak-overlap detection); `argparse` (CLI script); `knitr`, `rmarkdown` (vignettes); `testthat` (tests)

## Quick Start

```r
library(cascade)

config <- create_config(
  cell_types  = c("Mono", "DC", "NK", "B", "CD4_T", "CD8_T", "other_T"),
  chromosomes = paste0("chr", 1:22),
  file_patterns = list(
    # Per-cell-type ACAT significance
    eqtl_acat        = "data/eqtl/{CELL_TYPE}.acat.tsv.gz",
    caqtl_acat       = "data/caqtl/{CELL_TYPE}.acat.tsv.gz",
    # Per-cell-type, per-chromosome SuSiE fine-mapping
    eqtl_susie       = "data/eqtl/{CELL_TYPE}.{CHR}.susie.tsv.gz",
    caqtl_susie      = "data/caqtl/{CELL_TYPE}.{CHR}.susie.tsv.gz",
    # Mash LFSR matrices (one row per feature x variant, columns = cell types)
    eqtl_lfsr        = "data/eqtl.lfsr.tsv.gz",
    caqtl_lfsr       = "data/caqtl.lfsr.tsv.gz",
    # Meta-analysis results (Cochran's Q heterogeneity)
    eqtl_meta        = "data/eqtl.meta.tsv.gz",
    caqtl_meta       = "data/caqtl.meta.tsv.gz",
    # Peak coordinates and per-cell-type peak-gene links
    peak_bed         = "data/{CELL_TYPE}.peaks.bed",
    peak_gene_links  = "data/{CELL_TYPE}.peak_gene_links.tsv.gz"
  )
)

results <- run_cascade(config, output_dir = "results/")
```

See [`docs/SPECIFICATION.md`](docs/SPECIFICATION.md) §7 for the full input schema.

## Cell Type Hierarchy

CASCADE uses a configurable hierarchy to determine specificity categories. The package ships with `DEFAULT_CELL_HIERARCHY` for immune cells:

```r
DEFAULT_CELL_HIERARCHY
#> CellTypeHierarchy
#>   Lineages (2): myeloid, lymphoid
#>   L1 cell types (7): Mono, DC, NK, B, CD4_T, CD8_T, other_T
#>   T-cell subgroup: CD4_T, CD8_T, other_T
#>   Bulk: PBMC
#>   Categories (6): Cross-lineage shared | Likely shared but underpowered |
#>     Lineage-specific | T-cell-specific | Single cell-type | No significance
```

Define a custom hierarchy for any tissue (the category count adapts to the number of grouping levels you supply):

```r
brain_hierarchy <- create_cell_hierarchy(
  lineages = list(
    neuronal = c("Excitatory", "Inhibitory"),
    glial    = c("Astrocyte", "Oligodendrocyte", "Microglia", "OPC")
  ),
  bulk = "BulkBrain"
)
```

## Specificity Categories (default immune hierarchy)

| Category | Definition |
|---|---|
| Cross-lineage shared | Significant (FDR < 0.05) in ≥ 2 top-level lineages |
| Likely shared but underpowered | Apparent specificity demoted by LFSR gray-zone evidence (`0.05 ≤ LFSR < 0.5`) of hidden sharing |
| Lineage-specific | Significant in exactly one lineage, with high LFSR (≥ 0.5) elsewhere |
| T-cell-specific | Significant only within a defined subgroup (T-cell subgroup by default) |
| Single cell-type | Significant in exactly one L1 cell type |
| No significance | No significant effects detected |

Subgroup-specific categories are determined by the hierarchy. The default immune hierarchy yields one subgroup label (`T-cell-specific`); custom hierarchies may produce more.

## QTL Mechanisms

Each variant is classified along four features (chromatin accessibility peak overlap, caQTL status, peak-gene links, and eQTL status), yielding 25 mutually exclusive patterns that collapse into **8 mechanism categories** ordered by how directly the variant acts on its own regulatory sequence. The three cascade tiers are distinguished by where the causal chromatin accessibility peak sits relative to the variant.

| # | Mechanism | Definition |
|---:|---|---|
| 1 | Local Cascade | Variant overlaps a peak, shows caQTL effects for that overlapping peak, the peak has a peak-gene link, and the variant shows eQTL effects for that linked gene |
| 2 | Positional Cascade | Variant overlaps a peak with a peak-gene link and shows eQTL effects for the linked gene, but lacks detectable caQTL effects at the overlapping peak |
| 3 | Distal Cascade | Variant does not overlap any peak, but shows caQTL effects for a non-overlapping peak that has a peak-gene link to a gene for which the variant is also an eQTL |
| 4 | caQTL + eQTL (No Link) | Variant shows both caQTL and eQTL effects, but the eQTL gene differs from the linked gene or no peak-gene link exists |
| 5 | Only caQTL (With Link) | Variant shows caQTL effects and the affected peak has peak-gene links, but no eQTL is detected for the linked gene |
| 6 | Only caQTL (No Link) | Variant shows caQTL effects but the affected peak lacks peak-gene links |
| 7 | Only eQTL | Variant shows eQTL effects but no detectable caQTL effects |
| 8 | No molQTL | No significant molecular QTL effects in any cell type |

See `QTL_PATTERNS` and `QTL_MECHANISMS` for the full enumeration.

## Output

`run_cascade()` writes the following to `output_dir`:

| File | Contents |
|---|---|
| `gene_categorization.tsv.gz` | Per-gene specificity categories with heterogeneity-derived columns (top variant, max PIP, Cochran's Q) |
| `gene_categorization.summary.tsv` | Tabular summary of gene categorization |
| `peak_categorization.tsv.gz` | Per-peak specificity categories (same schema as gene output) |
| `peak_categorization.summary.tsv` | Tabular summary of peak categorization |
| `variant_categorization.tsv.gz` | Per-variant cross-cell-type aggregate: QTL mechanism + three specificity columns (combined / gene / peak) |
| `variant_categorization.qtl_mechanism.summary.tsv` | Counts of variants per QTL mechanism category |
| `variant_categorization.cell_type_specificity.summary.tsv` | Counts of variants per cell type specificity category |
| `variant_categorization_per_celltype_{CELL_TYPE}.tsv.gz` | Per-variant Stage 1 output, one file per cell type |

See [`docs/SPECIFICATION.md`](docs/SPECIFICATION.md) §8 for the full column schema.

## Citation

Kanai, M. et al. [Population-scale multiome immune cell atlas reveals complex disease drivers](https://doi.org/10.1101/2025.11.25.25340489). medRxiv (2025)

## License

MIT License

## Contact

Masahiro Kanai (<mkanai@broadinstitute.org>)
