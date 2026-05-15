# CASCADE Technical Specification

**Version**: 1.0.0

## 1. Overview

CASCADE (Comprehensive Analysis Suite for Cell-type specificity and Accessibility-Driven QTL Effects) is an R package that categorizes the cell type specificity of genes, chromatin accessibility peaks, and genetic variants from multi-cell-type QTL analyses. It addresses a fundamental confound in QTL studies: distinguishing genuine cell type-specific effects from apparent specificity caused by power differences across cell types (e.g., sample size, expression level, or chromatin accessibility variation).

The package implements two complementary analyses. **Feature-level analysis** (genes and peaks) assigns each gene or peak to one of N specificity categories based on ACAT significance across cell types and LFSR gray-zone evidence of hidden sharing. **Variant-level analysis** classifies each variant into one of 25 QTL patterns describing the mechanistic path from variant to phenotype (chromatin accessibility, gene expression, or cascade), then aggregates patterns across cell types to produce a final mechanism category and a secondary specificity assessment. The variant analysis itself runs in two stages (Section 6): per-cell-type pattern detection followed by cross-cell-type aggregation.

## 2. Cell Type Hierarchy

### 2.1 Hierarchy Structure

A `CellTypeHierarchy` object defines how L1 (primary) cell types are grouped for specificity categorization. It contains five components:

| Component | Type | Description |
|-----------|------|-------------|
| `lineages` | Named list of character vectors | 2+ top-level groups. Each L1 cell type belongs to exactly one lineage. |
| `subgroups` | Ordered list of named lists | 0+ nested grouping levels, ordered broadest to narrowest. Each subgroup must be a subset of a single lineage. |
| `bulk` | Character vector | Mixed/bulk cell types (e.g., PBMC). Bulk-only significance triggers "underpowered." |
| `other` | Character vector | Cell types outside all lineages. Participate in single-cell-type detection but not lineage counting. |
| `mapping_to_l1` | Named list | Maps lower-resolution types (L2, L3, ...) to their L1 parents. All categorization operates on L1. |

### 2.2 Dynamic Category Count

Total specificity categories = **4 fixed + N grouping levels**, where N = 1 (lineage level) + number of subgroup levels.

Fixed categories (always present):

1. Cross-lineage shared
2. Likely shared but underpowered
3. Single cell-type
4. No significance

Hierarchy-dependent categories (one per grouping level):

- `{Lineage label}-specific` (always present, broadest)
- `{Subgroup label}-specific` (one per subgroup level, narrowest last)

### 2.3 Default Immune Hierarchy

The shipped `DEFAULT_CELL_HIERARCHY` produces 6 categories (N=2 grouping levels):

```
Lineages:
  myeloid  = {Mono, DC}
  lymphoid = {NK, B, CD4_T, CD8_T, other_T}

Subgroups (1 level):
  T-cell = {CD4_T, CD8_T, other_T}

Bulk:  {PBMC}
Other: {other}

Categories: Cross-lineage shared | Likely shared but underpowered |
            Lineage-specific | T-cell-specific | Single cell-type | No significance
```

L2-to-L1 mapping includes 25 entries (e.g., CD14_Mono -> Mono, cDC1 -> DC, B_naive -> B, Treg -> CD4_T, ILC -> other).

### 2.4 Custom Hierarchies

Users define hierarchies via `create_cell_hierarchy()`. Example for brain tissue:

```r
brain_hierarchy <- create_cell_hierarchy(
  lineages = list(
    neuronal   = c("ExN", "InN"),
    glial      = c("Astro", "Oligo", "OPC", "Micro")
  ),
  subgroups = list(
    list(macroglial = c("Astro", "Oligo", "OPC"))
  ),
  bulk = "Cortex_bulk",
  column_prefix = "celltype"
)
# Produces 6 categories: Cross-lineage | Underpowered | Lineage-specific |
#   Macroglial-specific | Single cell-type | No significance
```

Constraints: L1 types must not overlap across lineages. Lower-level types may overlap freely (they are mapped to L1 before categorization). Dynamic per-feature hierarchies are not supported.

## 3. Cell Type Specificity Categories

### 3.1 Categorization Algorithm

Implemented in C++ (`categorization_core.cpp`). Input: significant cell types (from ACAT q < threshold), per-cell-type LFSR values, and the hierarchy's `CellTypeSets`. Steps:

```
1. Partition significant types into sig_core (non-bulk) and has_bulk flag.
2. If sig_core is empty and has_bulk -> "Likely shared but underpowered"
3. If sig_core is empty and no bulk  -> "No significance"
4. Count lineages with >= 1 type in sig_core.
   ("other" types are not in any lineage; they are skipped here.)
5. If 2+ lineages have signal -> "Cross-lineage shared"
   (No LFSR check. This is the broadest category; demotion is not applicable.)
6. If |sig_core| == 1 -> candidate "Single cell-type"
   LFSR check: ALL non-significant, non-bulk tested types (including "other").
   If gray zone found -> "Likely shared but underpowered"
   Else -> "Single cell-type"
7. From NARROWEST subgroup level to BROADEST:
   If all sig_core types fall within a single group at this level -> candidate "{group}-specific"
   LFSR check: sibling groups + other lineages (NOT same-group types, NOT "other" types).
   If gray zone found -> "Likely shared but underpowered"
   Else -> "{group}-specific"
   (First matching level wins; narrowest is checked first.)
8. If 1 lineage has signal (but no subgroup matched) -> candidate "Lineage-specific"
   LFSR check: types in OTHER lineages only (NOT same-lineage, NOT "other" types).
   If gray zone found -> "Likely shared but underpowered"
   Else -> "Lineage-specific"
```

### 3.2 LFSR Gray Zone Demotion

The gray zone is defined as `lfsr_sig_threshold <= LFSR < lfsr_null_threshold` (defaults: 0.05 <= LFSR < 0.5). A single cell type in the gray zone is sufficient to demote.

The scope of the LFSR check varies per candidate category (asymmetric by design):

| Candidate category | Types checked for gray zone | Types excluded | Rationale |
|---|---|---|---|
| Cross-lineage shared | None | All | Already broadest; no demotion possible |
| Lineage-specific | Other-lineage types | Same-lineage, "other" | Only cross-lineage evidence demotes |
| Subgroup-specific | Sibling subgroups + other lineages | Same-subgroup, "other" | Only evidence of broader sharing demotes |
| Single cell-type | All non-sig, non-bulk (incl. "other") | None | Most conservative; any evidence demotes |

This asymmetry is intentional. "Single cell-type" is the strongest specificity claim and requires the most conservative check. "Other" cell types (e.g., HSPC, ILC) have high gray-zone rates due to low power; including them in lineage/subgroup checks would collapse most categories into "underpowered."

### 3.3 LFSR Thresholds

| Parameter | Default | Operator | Meaning |
|-----------|---------|----------|---------|
| `lfsr_sig_threshold` | 0.05 | `<` | Significant effect |
| `lfsr_null_threshold` | 0.5 | `>=` | No evidence of effect |
| Gray zone | [0.05, 0.5) | `>= sig AND < null` | Ambiguous; potential hidden sharing |

## 4. QTL Mechanism Classification

### 4.1 Decision Variables

Each variant is classified per cell type using four boolean tests:

1. **Peak overlap**: Does the variant overlap a chromatin accessibility peak?
2. **caQTL**: Is the variant a significant caQTL? If so, for which peak (overlapping or non-overlapping)?
3. **Peak-gene link**: Does the caQTL peak have a significant link to a gene? (Overlapping link vs. non-overlapping link)
4. **eQTL**: Is the variant a significant eQTL? If so, for a linked gene or an unlinked gene?

### 4.2 The 25 QTL Patterns

The cascade mechanisms (Local / Positional / Distal) are distinguished by where
the chromatin accessibility bridge comes from. Within every mechanism group,
patterns are ordered so that overlap-peak gene links rank higher than
non-overlap-peak links, and peak-overlapping variants rank higher than
non-overlapping variants.

| Pat | Overlap | caQTL-ovl | caQTL-nonovl | Link-ovl | Link-nonovl | eQTL | eQTL-linked-ovl | eQTL-linked-nonovl | Interpretation | Mechanism |
|-----|---------|-----------|--------------|----------|-------------|------|------------------|---------------------|----------------|-----------|
| 1 | Y | Y | - | Y | - | Y | Y | - | Local Cascade | Local Cascade |
| 2 | Y | - | Y | Y | - | Y | Y | - | Positional Cascade (overlap link, non-overlap caQTL) | Positional Cascade |
| 3 | Y | N | N | Y | - | Y | Y | - | Positional Cascade (overlap link) | Positional Cascade |
| 4 | Y | - | Y | - | Y | Y | - | Y | Positional Cascade (non-overlap link) | Positional Cascade |
| 5 | N | - | Y | - | Y | Y | - | Y | Distal Cascade | Distal Cascade |
| 6 | Y | Y | - | Y | - | Y | N | - | caQTL + eQTL (overlap link, discordant) | caQTL + eQTL (No Link) |
| 7 | Y | - | Y | Y | - | Y | N | - | caQTL + eQTL (overlap link, non-overlap caQTL, discordant) | caQTL + eQTL (No Link) |
| 8 | Y | - | Y | - | Y | Y | - | N | caQTL + eQTL (non-overlap link, discordant) | caQTL + eQTL (No Link) |
| 9 | N | - | Y | - | Y | Y | - | N | caQTL + eQTL (no overlap, discordant) | caQTL + eQTL (No Link) |
| 10 | Y | Y | - | N | - | Y | - | - | caQTL + eQTL (no link) | caQTL + eQTL (No Link) |
| 11 | Y | - | Y | N | N | Y | - | - | caQTL + eQTL (non-overlap caQTL, no link) | caQTL + eQTL (No Link) |
| 12 | N | - | Y | - | N | Y | - | - | caQTL + eQTL (no overlap, no link) | caQTL + eQTL (No Link) |
| 13 | Y | Y | - | Y | - | N | - | - | Only caQTL (overlap link) | Only caQTL (With Link) |
| 14 | Y | - | Y | Y | - | N | - | - | Only caQTL (overlap link, non-overlap caQTL) | Only caQTL (With Link) |
| 15 | Y | - | Y | - | Y | N | - | - | Only caQTL (non-overlap link) | Only caQTL (With Link) |
| 16 | N | - | Y | - | Y | N | - | - | Only caQTL (no overlap, link) | Only caQTL (With Link) |
| 17 | Y | Y | - | N | - | N | - | - | Only caQTL (overlap, no link) | Only caQTL (No Link) |
| 18 | Y | - | Y | N | N | N | - | - | Only caQTL (non-overlap caQTL, no link) | Only caQTL (No Link) |
| 19 | N | - | Y | - | N | N | - | - | Only caQTL (no overlap, no link) | Only caQTL (No Link) |
| 20 | Y | N | N | Y | - | Y | N | - | Only eQTL (overlap link) | Only eQTL |
| 21 | Y | N | N | N | - | Y | - | - | Only eQTL (no link) | Only eQTL |
| 22 | N | N | N | - | - | Y | - | - | Only eQTL (no overlap) | Only eQTL |
| 23 | Y | N | N | Y | - | N | - | - | No molQTL (overlap link) | No molQTL |
| 24 | Y | N | N | N | - | N | - | - | No molQTL (no link) | No molQTL |
| 25 | N | N | N | - | - | N | - | - | No molQTL (no overlap) | No molQTL |

### 4.3 The 8 Mechanism Categories

| Index | Category | Patterns | Description |
|-------|----------|----------|-------------|
| 1 | Local Cascade | 1 | Variant overlaps peak, caQTL for that peak, link to gene, eQTL for linked gene |
| 2 | Positional Cascade | 2-4 | Variant overlaps peak; gene link and eQTL present but caQTL is absent for the linked peak or link runs through a different peak |
| 3 | Distal Cascade | 5 | Variant does not overlap any peak; caQTL for a distal peak, gene link, and eQTL all align |
| 4 | caQTL + eQTL (No Link) | 6-12 | Both caQTL and eQTL present but eQTL is discordant with the linked gene, or no peak-gene link exists |
| 5 | Only caQTL (With Link) | 13-16 | caQTL with a peak-gene link but no eQTL detected for that gene |
| 6 | Only caQTL (No Link) | 17-19 | caQTL without any peak-gene link |
| 7 | Only eQTL | 20-22 | eQTL without caQTL |
| 8 | No molQTL | 23-25 | No molecular QTL detected |

The mechanism hierarchy is ordered: lower index = more informative mechanism. Within each mechanism group, lower pattern number also reflects a more-informative configuration (overlap-peak link > non-overlap-peak link > no peak overlap).

## 5. Variant Heterogeneity

### 5.1 Heterogeneity Codes

| Code | Name | Description |
|------|------|-------------|
| a | shared_consistent | Same causal variant, consistent effect sizes across cell types |
| b | shared_heterogeneous | Same causal variant, significantly different effect magnitudes |
| c | shared_opposite | Same causal variant, opposite effect directions |
| d | distinct_variants | Different causal variants in different cell types |

### 5.2 Determination Methods

**Cochran's Q method** (default): Uses pre-computed meta-analysis results. Heterogeneity is significant when the Cochran's Q p-value falls below the configured threshold (default: 5e-8). If significant, effect direction determines code b (same direction) vs. c (opposite). If a feature has multiple variants, the method checks whether variants share credible sets.

**CS cluster method** (when configured): Uses pre-computed credible set cluster assignments. Variants in the same cluster across cell types are "shared" (codes a-c based on effect comparison). Variants in different clusters are "distinct" (code d). This method provides finer resolution than Cochran's Q by leveraging fine-mapping output directly.

### 5.3 Hierarchical Context for Variant Sharing

Whether a variant is "shared" depends on the gene's specificity category. A variant must have high PIP (>0.5) in at least 2 of the gene's significant cell types AND span the appropriate hierarchical level. The table below uses the default immune hierarchy for concreteness; custom hierarchies substitute their own lineage and subgroup labels (e.g., a brain hierarchy would replace "myeloid/lymphoid" with "neuronal/glial" and "T-cell-specific" with whatever subgroup the user defines).

| Gene specificity | Variant "shared" requirement |
|---|---|
| Cross-lineage shared | Variant must have high PIP in both myeloid AND lymphoid cell types |
| Lineage-specific | Variant must have high PIP in 2+ cell types within the same lineage |
| T-cell-specific (subgroup) | Variant must have high PIP in 2+ T-cell subtypes |
| Single cell-type / Underpowered | Any variant with high PIP in 2+ cell types is shared |

If no variant meets the sharing criteria, the feature is classified as `distinct_variants` (code d). If shared variants exist, Cochran's Q determines the heterogeneity sub-code:

1. **shared_consistent** (a): Q p-value above threshold (no significant heterogeneity)
2. **shared_heterogeneous** (b): Q p-value below threshold, same effect direction
3. **shared_opposite** (c): Q p-value below threshold, opposite effect directions (+- or -+)

When multiple shared variants exist, the most extreme pattern wins (opposite > heterogeneous > consistent). This hierarchical context ensures that a variant in B+NK cells is "shared" for a lymphoid-specific gene but NOT shared for a cross-lineage gene (since it does not span lineages).

**Underpowered leniency**: when the gene's specificity is "Likely shared but underpowered" (or the variant itself fails to span the gene's required hierarchical level), any variant present in 2+ cell types is treated as shared. This avoids over-classifying variants as `distinct_variants` purely because the underlying gene was demoted by the LFSR gray-zone check.

## 6. Two-Stage Variant Categorization

### 6.1 Stage 1: Per-Cell-Type QTL Pattern Detection

For each cell type independently:

1. Load SuSiE fine-mapping results (95% credible set variants passing PIP thresholds).
2. For each variant, determine: peak overlap (from BED intersection), caQTL status (SuSiE for peaks), peak-gene link presence (from pre-filtered link files), and eQTL status (SuSiE for genes).
3. Map the boolean combination to one of the 25 QTL patterns.
4. Record associated genes, peaks, and cascade peak-genes per variant.

**PIP filtering** is two-tiered:
- `pip_threshold` (default 0.5): max PIP across all cell types must exceed this (identity gate).
- `min_pip_threshold` (default 0.1): per-cell-type PIP must exceed this (association gate, intentionally lenient for underpowered cell types).

### 6.2 Stage 2: Cross-Cell-Type Aggregation

For each variant across all cell types:

1. **Best pattern**: `pmin(pattern_number)` across all cell types. Since patterns 23-25 (No molQTL) are the highest numbers, any cell type with patterns 1-22 automatically wins. For variants with no QTL evidence in any cell type, the best pattern is 25.
2. **Mechanism assignment**: `QTL_PATTERNS[[best_pattern]]$mechanism` maps pattern to one of 8 categories.
3. **Cell type sets**:
   - `significant_cts`: all cell types with any QTL evidence (patterns 1-22).
   - `gene_affected_cell_types`: cell types with patterns in {1-12, 20-22} (has eQTL component).
   - `peak_affected_cell_types`: cell types with patterns in {1-2, 4-19} (has caQTL component; pattern 3 is excluded because it lacks a caQTL).
4. **Best cell types**: cell types where the pattern equals the best pattern.
5. **Three-way specificity**: Each variant receives three specificity assessments:
   - `cell_type_specificity`: based on cell types sharing the best mechanism category AND having QTL evidence. LFSR disabled.
   - `gene_cell_type_specificity`: based on `gene_affected_cell_types` only. LFSR enabled (eQTL LFSR).
   - `peak_cell_type_specificity`: based on `peak_affected_cell_types` only. LFSR enabled (caQTL LFSR).
6. **No synthetic mechanisms**: If different cell types have different mechanism categories, only the best (lowest pattern number) is reported.

## 7. Data Requirements

### 7.1 Input File Types

The `file_patterns` slot of the config maps each input to a path or path template. Templates may include `{CELL_TYPE}` and `{CHR}` placeholders that the loader substitutes at runtime.

| `file_patterns` key | Contents | Required columns (defaults) |
|---------------------|----------|-----------------------------|
| `eqtl_acat` | Gene-level ACAT significance per cell type | phenotype_id, ACAT_q |
| `caqtl_acat` | Peak-level ACAT significance per cell type | phenotype_id, qval |
| `eqtl_susie` | Gene-level SuSiE fine-mapping per cell type | rsid, region, prob, chromosome, cs, beta, se |
| `caqtl_susie` | Peak-level SuSiE fine-mapping per cell type | rsid, region, prob, chromosome, cs, beta, se |
| `peak_gene_links` | Peak-to-gene links (pre-filtered) | peak_id, gene_id |
| `peak_bed` | BED of peak coordinates per cell type (used for variant overlap) | chrom, start, end, peak_id |
| `eqtl_lfsr` / `caqtl_lfsr` | Mash LFSR values (wide format, one column per cell type) | id (with separator) |
| `eqtl_meta` / `caqtl_meta` | Meta-analysis results (shared schema) | variant, phenotype, meta_nlog10p_het, direction, max_pip, max_chisq |
| `eqtl_mashr` / `caqtl_mashr` | mashr model RDS files (optional; used for L1/L2 sharing analysis) | — |
| `cs_clusters` | Credible set cluster assignments (optional) | ID, QTL, cell_type, trait, cs, cluster |
| `cs_cluster_variants` | Cluster-level variant details (optional) | cluster, variant_id, features |

**Naming note:** the `column_mapping` slot uses unprefixed keys (`peak_gene`, `lfsr`, `meta`) because the column schema is the same for the eQTL and caQTL variants of those files. Only `file_patterns` distinguishes them with `eqtl_`/`caqtl_` prefixes (since each modality has a separate file path).

### 7.2 Column Mapping

All input column names are configurable via `DEFAULT_COLUMN_MAPPING` or user overrides in `create_config(column_mapping = ...)`. The mapping is a nested list keyed by file type, then `internal_name = "input_column"`. User overrides are merged per-file-type with defaults.

**DEFAULT_COLUMN_MAPPING schema** (9 file types):

| File type | Internal name | Default input column | Notes |
|-----------|--------------|---------------------|-------|
| `eqtl_acat` | feature_id, q_value | phenotype_id, ACAT_q | |
| `caqtl_acat` | feature_id, q_value | phenotype_id, qval | |
| `eqtl_susie` | variant_id, feature_id, pip, chromosome, cs_id, beta, se | rsid, region, prob, chromosome, cs, beta, se | beta/se optional (NULL to skip) |
| `caqtl_susie` | variant_id, feature_id, pip, chromosome, cs_id, beta, se | rsid, region, prob, chromosome, cs, beta, se | beta/se optional |
| `peak_gene` | peak_id, gene_id | peak_id, gene_id | |
| `lfsr` | id_column, id_separator | id, : | ID split into feature_id + variant_id inside `load_and_parse_lfsr_file()` using the configured separator |
| `meta` | variant_id, feature_id, cochran_q_nlog10p, direction, max_pip, max_chisq | variant, phenotype, meta_nlog10p_het, direction, max_pip, max_chisq | cochran_q_pval derived via `10^(-nlog10p)` |
| `cs_clusters` | variant_id, qtl_type, cell_type, feature_id, cs_id, cluster_id | ID, QTL, cell_type, trait, cs, cluster | |
| `cs_cluster_variants` | cluster_id, variant_id, feature_ids | cluster, variant_id, features | feature_ids is comma-separated |

**Column mapping utilities:**

- `rename_columns(dt, mapping)`: 1:1 rename from input to internal names; NULL-valued entries are optional and skipped; errors on missing required columns.
- `load_and_parse_lfsr_file(path, mapping)`: Loads an LFSR file and splits the composite ID column (e.g., `gene:variant`) into `feature_id` and `variant_id`.
- `detect_ct_columns(dt, hierarchy)`: Uses the cell type hierarchy to identify which columns in LFSR/ACAT data correspond to cell types.

### 7.3 Configuration

Configuration is provided as an R list (via `create_config()`) or a JSON file. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pip_threshold` | 0.5 | Max PIP across cell types for variant inclusion |
| `min_pip_threshold` | 0.1 | Per-cell-type PIP threshold |
| `acat_fdr_threshold` | 0.05 | ACAT FDR for variant data loading |
| `lfsr_sig_threshold` | 0.05 | LFSR significance threshold |
| `lfsr_null_threshold` | 0.5 | LFSR null threshold |
| `cochran_q_threshold` | 5e-8 | Cochran's Q p-value threshold for heterogeneity |
| `n_cores` | auto | Parallelization cores |
| `chromosomes` | _required_ | Chromosomes to analyze (no default; loader aborts if unset) |

The `hierarchy` field accepts a `CellTypeHierarchy` object or a nested list (auto-resolved). If absent, `DEFAULT_CELL_HIERARCHY` is used.

## 8. Output Schema

### 8.1 Gene/Peak Categorization Output

One row per feature. Columns produced by `categorize_features()` (via its underlying `categorize_features_from_acat` C++ kernel):

| Column | Type | Description |
|--------|------|-------------|
| `feature_id` | character | Feature identifier (gene or peak ID) |
| `cell_type_specificity` | character | Assigned specificity category (Section 3) |
| `cell_type_specificity_pattern` | character | Pattern index (1..N categories) used internally for ordering |
| `significant_cts` | character | Comma-separated significant L1 cell types (deterministic order) |
| `tested_cts` | character | Comma-separated tested L1 cell types (includes bulk) |
| `n_significant_cts` | integer | Number of significant L1 cell types |
| `n_tested_cts` | integer | Number of tested L1 cell types |

When variant analysis runs (heterogeneity / linkage), the following columns are appended:

| Column | Type | Description |
|--------|------|-------------|
| `variant_heterogeneity` | character | Heterogeneity code (Section 5.1) |
| `variant_heterogeneity_pattern` | character | Pattern index used internally |
| `n_variants` | integer | Number of variants considered |
| `cs_cts` | character | Cell types contributing credible sets |
| `n_cs_cts` | integer | Number of CS-contributing cell types |
| `top_variant` | character | Lead variant ID |
| `max_pip` | numeric | Maximum PIP across cell types for the lead variant |
| `max_chisq` | numeric | Maximum chi-square across cell types |
| `hierarchical_variant_pattern` | character | Pattern reflecting the hierarchical sharing context (Section 5.3) |

### 8.2 Variant Cross-Cell-Type Output

One row per variant. 19 columns (matching `final_results` in `variant_stage2_aggregation.R`):

| Column | Type | Description |
|--------|------|-------------|
| `variant_id` | character | Variant identifier (rsID) |
| `qtl_mechanism_category` | character | One of 8 mechanism categories (Section 4.3) |
| `cell_type_specificity` | character | Combined specificity based on mechanism-sharing cell types (LFSR disabled) |
| `gene_cell_type_specificity` | character | Gene-specific specificity (eQTL cell types only, with LFSR) |
| `peak_cell_type_specificity` | character | Peak-specific specificity (caQTL cell types only, with LFSR) |
| `qtl_pattern_number` | integer | Best pattern number (1-25) |
| `qtl_pattern` | character | Human-readable pattern interpretation |
| `best_cell_types` | character | Cell types with the best (lowest) pattern |
| `significant_cts` | character | All cell types with any QTL evidence (patterns 1-22) |
| `gene_affected_cell_types` | character | Cell types with eQTL component |
| `peak_affected_cell_types` | character | Cell types with caQTL component |
| `associated_genes` | character | Comma-separated eQTL target genes |
| `associated_peaks` | character | Comma-separated caQTL target peaks |
| `cascade_peak_genes` | character | Genes linked via cascade (peak -> gene) |
| `peak_overlap` | logical | Whether variant overlaps any peak (from best cell type) |
| `caqtl` | character | caQTL status description |
| `peak_gene_link` | character | Peak-gene link status |
| `eqtl` | character | eQTL status description |
| `link_mechanism` | character | Link mechanism between peak and gene |

### 8.3 Variant Per-Cell-Type Output

One file per cell type. Contains the Stage 1 results from `categorize_qtl_patterns_cpp`:
`variant_id`, `qtl_pattern_number`, `qtl_mechanism_category`, `associated_genes`, `associated_peaks`, `cascade_peak_genes`, `peak_overlap`, `caqtl`, `peak_gene_link`, `eqtl`, `link_mechanism`.

### 8.4 Field Ordering Guarantees

R-side comma-separated fields (`best_cell_types`, `significant_cts`, `gene_affected_cell_types`, `peak_affected_cell_types`) are deterministic, ordered by cell type column position. `link_mechanism` is also deterministic (built from sorted unique values). C++-side fields (`associated_genes`, `associated_peaks`, `cascade_peak_genes`) use unordered sets and have non-deterministic ordering across runs.
