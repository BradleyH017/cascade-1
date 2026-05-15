#' @title Cell Type Hierarchy
#' @description Functions for defining and working with cell type hierarchies
#'   used in specificity categorization.

#' Create a Cell Type Hierarchy
#'
#' Defines the cell type hierarchy for specificity categorization. Each grouping
#' level adds one specificity category. The package ships with
#' \code{DEFAULT_CELL_HIERARCHY} for immune cells (2 grouping levels -> 6
#' categories), but users can define their own system for any tissue or organism.
#'
#' The categorization logic uses the hierarchy as follows:
#' \itemize{
#'   \item "Cross-lineage shared": significant in 2+ top-level lineage groups
#'   \item "Likely shared but underpowered": LFSR gray zone evidence of hidden sharing
#'   \item "{Lineage}-specific": significant in exactly 1 lineage group
#'   \item "{Subgroup}-specific": significant in 1 subgroup (one category per subgroup level)
#'   \item "Single cell-type": significant in exactly 1 L1 cell type
#'   \item "No significance": nothing significant
#' }
#'
#' Total categories = 4 fixed + N grouping levels (1 lineage level + subgroup levels).
#'
#' @param lineages Named list with 2+ entries. Each entry is a character vector
#'   of L1 cell type names in that top-level group. Names become lineage labels.
#' @param subgroups Ordered list of sub-grouping levels (optional). Each element
#'   is a named list of groups at that depth. Groups must be subsets of a single
#'   lineage. Ordered from broadest to narrowest. Use \code{attr(level, "label")}
#'   to set a custom category label for a level with multiple groups.
#' @param bulk Character vector of bulk/mixed cell type names (e.g., "PBMC").
#'   Features significant only in bulk are categorized as "Likely shared but
#'   underpowered". NULL if none.
#' @param other Character vector of cell types excluded from lineage grouping
#'   (e.g., "other" for unclassified types). These types still participate in
#'   single-cell-type detection and are valid mapping targets. NULL if none.
#' @param mapping_to_l1 Named list mapping lower-level cell types to their L1
#'   parents. Supports arbitrary QTL resolution depth (L2, L3, etc.) -- all map
#'   directly to L1. Targets must be L1 lineage types or "other" types. Types
#'   not in this mapping are assumed to already be L1.
#' @param column_prefix Prefix for cell type column names in data files.
#'   Default: "predicted.celltype". Columns will be "{prefix}.l1.{name}" for L1.
#' @return A CellTypeHierarchy object (S3 class)
#' @export
create_cell_hierarchy <- function(lineages,
                                  subgroups = list(),
                                  bulk = NULL,
                                  other = NULL,
                                  mapping_to_l1 = list(),
                                  column_prefix = "predicted.celltype") {
  # --- Validation ---
  if (!is.list(lineages) || length(lineages) < 2 || is.null(names(lineages))) {
    stop("lineages must be a named list with 2+ entries")
  }

  all_l1 <- unlist(lineages, use.names = FALSE)
  if (anyDuplicated(all_l1)) {
    stop(
      "L1 cell types must not overlap across lineages: ",
      paste(all_l1[duplicated(all_l1)], collapse = ", ")
    )
  }

  # Validate subgroup levels
  for (i in seq_along(subgroups)) {
    level <- subgroups[[i]]
    if (!is.list(level) || is.null(names(level))) {
      stop("subgroups[[", i, "]] must be a named list")
    }
    for (nm in names(level)) {
      bad <- setdiff(level[[nm]], all_l1)
      if (length(bad) > 0) {
        stop(
          "Subgroup '", nm, "' contains cell types not in any lineage: ",
          paste(bad, collapse = ", ")
        )
      }
      # Check subgroup is within a single lineage
      parent_lineages <- vapply(lineages, function(lin) {
        all(level[[nm]] %in% lin)
      }, logical(1))
      if (!any(parent_lineages)) {
        stop("Subgroup '", nm, "' spans multiple lineages")
      }
    }
  }

  # Validate mapping targets
  valid_targets <- c(all_l1, other)
  for (target in unique(unlist(mapping_to_l1))) {
    if (!target %in% valid_targets) {
      stop(
        "mapping_to_l1 target '", target,
        "' is not an L1 cell type or 'other' type"
      )
    }
  }

  # Validate bulk/other are disjoint from lineages
  if (!is.null(bulk) && any(bulk %in% all_l1)) {
    stop("bulk types must not overlap with lineage L1 types")
  }
  if (!is.null(other) && any(other %in% all_l1)) {
    stop("other types must not overlap with lineage L1 types")
  }

  # --- Auto-generate category labels ---
  n_subgroup_levels <- length(subgroups)
  n_grouping_levels <- 1 + n_subgroup_levels

  level_labels <- "Lineage"
  subgroup_counter <- 0L
  for (level in subgroups) {
    if (!is.null(attr(level, "label"))) {
      level_labels <- c(level_labels, attr(level, "label"))
    } else if (length(level) == 1) {
      level_labels <- c(level_labels, names(level))
    } else {
      subgroup_counter <- subgroup_counter + 1L
      label <- if (subgroup_counter == 1L) "Subgroup" else paste0("Subgroup-L", subgroup_counter)
      level_labels <- c(level_labels, label)
    }
  }

  category_labels <- c(
    "Cross-lineage shared",
    "Likely shared but underpowered",
    paste0(level_labels, "-specific"),
    "Single cell-type",
    "No significance"
  )

  # --- Build full column names ---
  make_col <- function(level_tag, name) {
    paste0(column_prefix, ".", level_tag, ".", name)
  }
  l1_columns <- setNames(
    vapply(all_l1, function(n) make_col("l1", n), character(1)),
    all_l1
  )
  bulk_columns <- if (!is.null(bulk)) {
    setNames(vapply(bulk, function(n) make_col("l1", n), character(1)), bulk)
  } else {
    character(0)
  }
  other_columns <- if (!is.null(other)) {
    setNames(vapply(other, function(n) make_col("l1", n), character(1)), other)
  } else {
    character(0)
  }

  structure(list(
    lineages = lineages,
    subgroups = subgroups,
    bulk = bulk,
    other = other,
    mapping_to_l1 = mapping_to_l1,
    column_prefix = column_prefix,
    category_labels = category_labels,
    l1_columns = l1_columns,
    bulk_columns = bulk_columns,
    other_columns = other_columns,
    l1_cell_types = all_l1,
    n_grouping_levels = n_grouping_levels
  ), class = "CellTypeHierarchy")
}

#' Print method for CellTypeHierarchy
#' @param x A CellTypeHierarchy object
#' @param ... Additional arguments (ignored)
#' @export
print.CellTypeHierarchy <- function(x, ...) {
  cat("CellTypeHierarchy\n")
  cat("  Lineages (", length(x$lineages), "): ",
    paste(names(x$lineages), collapse = ", "), "\n",
    sep = ""
  )
  cat("  L1 cell types (", length(x$l1_cell_types), "): ",
    paste(x$l1_cell_types, collapse = ", "), "\n",
    sep = ""
  )
  if (length(x$subgroups) > 0) {
    cat("  Subgroup levels: ", length(x$subgroups), "\n", sep = "")
    for (i in seq_along(x$subgroups)) {
      cat("    Level ", i, ": ",
        paste(names(x$subgroups[[i]]), collapse = ", "), "\n",
        sep = ""
      )
    }
  }
  if (!is.null(x$bulk)) cat("  Bulk: ", paste(x$bulk, collapse = ", "), "\n", sep = "")
  if (!is.null(x$other)) cat("  Other: ", paste(x$other, collapse = ", "), "\n", sep = "")
  if (length(x$mapping_to_l1) > 0) {
    cat("  Mapping to L1: ", length(x$mapping_to_l1), " entries\n", sep = "")
  }
  cat("  Categories (", length(x$category_labels), "): ",
    paste(x$category_labels, collapse = " | "), "\n",
    sep = ""
  )
  invisible(x)
}

#' Default Cell Type Hierarchy (Immune)
#'
#' The default hierarchy for immune cell QTL analysis, with myeloid/lymphoid
#' lineages and T-cell subgroup. Produces 6 specificity categories.
#' @export
DEFAULT_CELL_HIERARCHY <- create_cell_hierarchy(
  lineages = list(
    myeloid  = c("Mono", "DC"),
    lymphoid = c("NK", "B", "CD4_T", "CD8_T", "other_T")
  ),
  subgroups = list(
    list("T-cell" = c("CD4_T", "CD8_T", "other_T"))
  ),
  bulk = "PBMC",
  other = "other",
  mapping_to_l1 = list(
    CD14_Mono = "Mono", CD16_Mono = "Mono",
    cDC1 = "DC", cDC2 = "DC", pDC = "DC",
    B_intermediate = "B", B_memory = "B", B_naive = "B", Plasmablast = "B",
    CD4_CTL = "CD4_T", CD4_Naive = "CD4_T", CD4_TCM = "CD4_T",
    CD4_TEM = "CD4_T", Treg = "CD4_T",
    CD8_Naive = "CD8_T", CD8_TEM = "CD8_T",
    NK = "NK", NK_CD56bright = "NK", NK_Proliferating = "NK",
    MAIT = "other_T", dnT = "other_T", gdT = "other_T",
    ILC = "other", HSPC = "other", Platelet = "other"
  ),
  column_prefix = "predicted.celltype"
)

#' Resolve Hierarchy in Config
#'
#' Normalizes the hierarchy field in a config object. If the hierarchy is
#' already a CellTypeHierarchy, returns as-is. If it's a plain list (e.g.,
#' from JSON), constructs a CellTypeHierarchy. If missing, uses the default.
#'
#' @param config Configuration list
#' @return Config with resolved hierarchy
#' @keywords internal
resolve_hierarchy <- function(config) {
  if (!inherits(config$hierarchy, "CellTypeHierarchy")) {
    if (is.list(config$hierarchy)) {
      config$hierarchy <- do.call(create_cell_hierarchy, config$hierarchy)
    } else {
      config$hierarchy <- DEFAULT_CELL_HIERARCHY
    }
  }
  # Backfill cell_types if not already set
  if (is.null(config$cell_types)) {
    config$cell_types <- config$hierarchy$l1_cell_types
  }
  config
}

#' Extract C++ Interface Parameters from Hierarchy
#'
#' Bridge function: converts a CellTypeHierarchy into the parameter format
#' expected by the C++ categorization functions (lineage_groups, subgroup_levels,
#' bulk_cts, other_cts, l2_to_l1_mapping, specificity_categories).
#'
#' @param hierarchy A CellTypeHierarchy object
#' @return Named list with C++ parameter values
#' @keywords internal
hierarchy_to_cpp_params <- function(hierarchy) {
  # Lineage groups: list of character vectors (full column names)
  lineage_groups <- lapply(hierarchy$lineages, function(types) {
    unname(hierarchy$l1_columns[types])
  })

  # Subgroup levels: list of lists of character vectors
  subgroup_levels <- lapply(hierarchy$subgroups, function(level) {
    lapply(level, function(types) {
      unname(hierarchy$l1_columns[types])
    })
  })

  # Bulk and other: character vectors of full column names
  bulk_cts <- unname(hierarchy$bulk_columns)
  other_cts <- unname(hierarchy$other_columns)

  list(
    lineage_groups = lineage_groups,
    subgroup_levels = subgroup_levels,
    bulk_cts = bulk_cts,
    other_cts = other_cts,
    l2_to_l1_mapping = get_mapping_columns(hierarchy),
    specificity_categories = hierarchy$category_labels
  )
}

#' Get Lineage Column Names
#'
#' Returns full column names for lineage groups from a hierarchy.
#' @param hierarchy A CellTypeHierarchy object
#' @return List of character vectors (one per lineage)
#' @keywords internal
get_lineage_columns <- function(hierarchy) {
  lapply(hierarchy$lineages, function(types) {
    hierarchy$l1_columns[types]
  })
}

#' Get Subgroup Level Column Names
#'
#' Returns full column names for each subgroup level.
#' @param hierarchy A CellTypeHierarchy object
#' @return List of lists of character vectors
#' @keywords internal
get_subgroup_level_columns <- function(hierarchy) {
  lapply(hierarchy$subgroups, function(level) {
    lapply(level, function(types) {
      hierarchy$l1_columns[types]
    })
  })
}

#' Get Mapping to L1 with Full Column Names
#'
#' Returns the mapping_to_l1 with full column name keys and values.
#' @param hierarchy A CellTypeHierarchy object
#' @return Named list mapping full lower-level column names to full L1 column names
#' @keywords internal
get_mapping_columns <- function(hierarchy) {
  if (length(hierarchy$mapping_to_l1) == 0) {
    return(list())
  }
  # Build full column names for mapping entries
  prefix <- hierarchy$column_prefix
  result <- list()
  for (lower_name in names(hierarchy$mapping_to_l1)) {
    l1_name <- hierarchy$mapping_to_l1[[lower_name]]
    lower_col <- paste0(prefix, ".l2.", lower_name)
    l1_col <- if (l1_name %in% names(hierarchy$l1_columns)) {
      hierarchy$l1_columns[[l1_name]]
    } else if (l1_name %in% names(hierarchy$other_columns)) {
      hierarchy$other_columns[[l1_name]]
    } else {
      paste0(prefix, ".l1.", l1_name)
    }
    result[[lower_col]] <- l1_col
  }
  result
}
