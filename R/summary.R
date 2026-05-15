#' @title Categorization Summary Functions
#' @description Generate per-column summary tables for gene, peak, and variant categorization results

#' Summarize a single categorical column
#'
#' @param data Data frame
#' @param column_name Name of column to summarize
#' @param header Optional header label inserted as the first row
#' @return Summary data frame with columns category, count, percentage
#' @keywords internal
summarize_column <- function(data, column_name, header = NULL) {
  if (!column_name %in% names(data)) {
    return(data.frame(
      category = character(0),
      count = integer(0),
      percentage = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  value_counts <- table(data[[column_name]], useNA = "ifany")
  total <- sum(value_counts)

  summary_df <- data.frame(
    category = names(value_counts),
    count = as.integer(value_counts),
    percentage = round(100 * as.numeric(value_counts) / total, 2),
    stringsAsFactors = FALSE
  )

  summary_df$category[is.na(summary_df$category)] <- "Not categorized"
  summary_df <- summary_df[order(summary_df$count, decreasing = TRUE), ]

  if (!is.null(header)) {
    header_row <- data.frame(
      category = header,
      count = total,
      percentage = 100.00,
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(header_row, summary_df)
  }

  total_row <- data.frame(
    category = "TOTAL",
    count = total,
    percentage = 100.00,
    stringsAsFactors = FALSE
  )
  summary_df <- rbind(summary_df, total_row)

  return(summary_df)
}

#' Save categorization summary
#'
#' @param summary_data Summary data frame
#' @param output_file Output file path
#' @import data.table
#' @keywords internal
save_categorization_summary <- function(summary_data, output_file) {
  data.table::fwrite(summary_data, output_file,
    sep = "\t", quote = FALSE, na = "NA"
  )
  cli::cli_alert_success("Saved summary to: {.file {output_file}}")
}
