#' @title Scientific Definitions
#' @description QTL mechanism categories, pattern definitions, variant heterogeneity
#'   codes, and other scientific constants used in categorization.

#' QTL Mechanism Categories
#'
#' Eight main categories for variant QTL mechanisms
#' @export
QTL_MECHANISMS <- c(
  "Local Cascade", # 1
  "Positional Cascade", # 2
  "Distal Cascade", # 3
  "caQTL + eQTL (No Link)", # 4 (includes discordant targets)
  "Only caQTL (With Link)", # 5
  "Only caQTL (No Link)", # 6
  "Only eQTL", # 7
  "No molQTL" # 8
)

#' QTL Pattern Details
#'
#' Maps 25 QTL patterns to their interpretations and mechanism categories
#' Each pattern has: interpretation (detailed description) and mechanism (category index)
#' @export
QTL_PATTERNS <- list(
  # Local Cascade (1)
  list(interpretation = "Local Cascade", mechanism = 1),

  # Positional Cascade (2-4): overlap link prioritized over non-overlap link
  list(interpretation = "Positional Cascade (overlap link, non-overlap caQTL)", mechanism = 2),
  list(interpretation = "Positional Cascade (overlap link)", mechanism = 2),
  list(interpretation = "Positional Cascade (non-overlap link)", mechanism = 2),

  # Distal Cascade (5)
  list(interpretation = "Distal Cascade", mechanism = 3),

  # caQTL + eQTL (No Link) (6-12): discordant-with-link first (by link type), then no-link (by caQTL/overlap type)
  list(interpretation = "caQTL + eQTL (overlap link, discordant)", mechanism = 4),
  list(interpretation = "caQTL + eQTL (overlap link, non-overlap caQTL, discordant)", mechanism = 4),
  list(interpretation = "caQTL + eQTL (non-overlap link, discordant)", mechanism = 4),
  list(interpretation = "caQTL + eQTL (no overlap, discordant)", mechanism = 4),
  list(interpretation = "caQTL + eQTL (no link)", mechanism = 4),
  list(interpretation = "caQTL + eQTL (non-overlap caQTL, no link)", mechanism = 4),
  list(interpretation = "caQTL + eQTL (no overlap, no link)", mechanism = 4),

  # Only caQTL (With Link) (13-16)
  list(interpretation = "Only caQTL (overlap link)", mechanism = 5),
  list(interpretation = "Only caQTL (overlap link, non-overlap caQTL)", mechanism = 5),
  list(interpretation = "Only caQTL (non-overlap link)", mechanism = 5),
  list(interpretation = "Only caQTL (no overlap, link)", mechanism = 5),

  # Only caQTL (No Link) (17-19)
  list(interpretation = "Only caQTL (overlap, no link)", mechanism = 6),
  list(interpretation = "Only caQTL (non-overlap caQTL, no link)", mechanism = 6),
  list(interpretation = "Only caQTL (no overlap, no link)", mechanism = 6),

  # Only eQTL (20-22)
  list(interpretation = "Only eQTL (overlap link)", mechanism = 7),
  list(interpretation = "Only eQTL (no link)", mechanism = 7),
  list(interpretation = "Only eQTL (no overlap)", mechanism = 7),

  # No molQTL (23-25)
  list(interpretation = "No molQTL (overlap link)", mechanism = 8),
  list(interpretation = "No molQTL (no link)", mechanism = 8),
  list(interpretation = "No molQTL (no overlap)", mechanism = 8)
)

#' Get QTL Pattern Interpretation
#'
#' @param pattern_num Numeric pattern number (1-25)
#' @return Character string with the pattern interpretation
get_pattern_interpretation <- function(pattern_num) {
  if (!is.numeric(pattern_num) || pattern_num < 1 || pattern_num > length(QTL_PATTERNS)) {
    stop(sprintf(
      "Invalid pattern number: %s. Must be between 1 and %d.",
      pattern_num, length(QTL_PATTERNS)
    ))
  }
  return(QTL_PATTERNS[[pattern_num]]$interpretation)
}

#' Variant Heterogeneity Categories
#'
#' Categories for variant heterogeneity in multi-cell-type features.
#' Maps category names to letter codes (a-d)
#' @export
VARIANT_HETEROGENEITY <- c(
  shared_consistent = "a", # Same variant, consistent effects
  shared_heterogeneous = "b", # Same variant, different magnitudes
  shared_opposite = "c", # Same variant, opposite directions
  distinct_variants = "d" # Different variants in different cells
)

#' LFSR Significance Threshold
LFSR_SIG_THRESHOLD <- 0.05

#' LFSR Null Hypothesis Threshold
LFSR_NULL_THRESHOLD <- 0.5

#' caQTL Status Descriptions
CAQTL_STATUS <- list(
  OVERLAPPING = "For overlapping peak",
  NON_OVERLAPPING = "For non-overlapping peak",
  NO_CAQTL = "No detected caQTL"
)

#' eQTL Status Descriptions
EQTL_STATUS <- list(
  LINKED = "For linked gene",
  NON_LINKED = "For non-linked gene",
  EQTL = "eQTL",
  NO_EQTL = "No detected eQTL"
)

#' Default Column Mappings
#'
#' Maps internal column names to input column names for each file type.
#' Entries with NULL values are optional and skipped during rename.
#' Users can override individual entries via create_config(column_mapping = list(...)).
DEFAULT_COLUMN_MAPPING <- list(
  eqtl_acat = list(
    feature_id = "phenotype_id",
    q_value = "ACAT_q"
  ),
  caqtl_acat = list(
    feature_id = "phenotype_id",
    q_value = "qval"
  ),
  eqtl_susie = list(
    variant_id = "rsid",
    feature_id = "region",
    pip = "prob",
    chromosome = "chromosome",
    cs_id = "cs",
    beta = "beta",
    se = "se"
  ),
  caqtl_susie = list(
    variant_id = "rsid",
    feature_id = "region",
    pip = "prob",
    chromosome = "chromosome",
    cs_id = "cs",
    beta = "beta",
    se = "se"
  ),
  peak_gene = list(
    peak_id = "peak_id",
    gene_id = "gene_id"
  ),
  lfsr = list(
    id_column = "id",
    id_separator = ":"
  ),
  meta = list(
    variant_id = "variant",
    feature_id = "phenotype",
    cochran_q_nlog10p = "meta_nlog10p_het",
    direction = "direction",
    max_pip = "max_pip",
    max_chisq = "max_chisq"
  ),
  cs_clusters = list(
    variant_id = "ID",
    qtl_type = "QTL",
    cell_type = "cell_type",
    feature_id = "trait",
    cs_id = "cs",
    cluster_id = "cluster"
  ),
  cs_cluster_variants = list(
    cluster_id = "cluster",
    variant_id = "variant_id",
    feature_ids = "features"
  )
)

#' Check if a cell type is a mapped (lower-level) type
#'
#' A cell type is considered "L2" (or lower) if it appears as a key
#' in the hierarchy's mapping_to_l1, meaning it maps to an L1 parent.
#'
#' @param celltype Character scalar or vector of cell type column names
#' @param hierarchy A CellTypeHierarchy object (required)
#' @return Logical vector
is_l2_celltype <- function(celltype, hierarchy) {
  mapped_cols <- names(get_mapping_columns(hierarchy))
  celltype %in% mapped_cols
}

#' Filter to only L1 cell types
#'
#' @param celltypes Character vector of cell type column names
#' @param hierarchy A CellTypeHierarchy object (required)
#' @return Character vector of L1-only cell types
filter_l1_celltypes <- function(celltypes, hierarchy) {
  celltypes[!is_l2_celltype(celltypes, hierarchy)]
}
