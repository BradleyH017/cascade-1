# test-specificity-categorization.R
#
# Tests every branch of the 6-category specificity decision tree using
# both the oracle (helper-oracle.R) and the C++ implementation
# (categorize_cell_specificity via cascade:::).
#
# 18 test cases covering:
#   - Cross-lineage shared (cases 1, 2, 13)
#   - Likely shared but underpowered (cases 3, 4, 8, 10, 14)
#   - Lineage-specific (cases 6, 7, 18)
#   - T-cell-specific (cases 9, 11)
#   - Single cell-type (cases 5, 15, 16, 17)
#   - No significance (case 12)

# Helper: call the C++ function with oracle-compatible arguments
call_cpp <- function(sig_l1_cts, lfsr_values = NULL) {
  lfsr_vals <- NULL
  lfsr_nms <- NULL
  if (!is.null(lfsr_values)) {
    lfsr_vals <- unname(lfsr_values)
    lfsr_nms <- names(lfsr_values)
  }
  cascade:::categorize_cell_specificity(
    significant_cts = sig_l1_cts,
    specificity_categories = cascade::DEFAULT_CELL_HIERARCHY$category_labels,
    lfsr_values = lfsr_vals,
    lfsr_names = lfsr_nms,
    lineage_groups = list(ORACLE_MYELOID, ORACLE_LYMPHOID),
    subgroup_levels = list(list(ORACLE_TCELL)),
    bulk_cts = ORACLE_PBMC,
    other_cts = "predicted.celltype.l1.other",
    lfsr_sig_threshold = 0.05,
    lfsr_null_threshold = 0.5
  )
}


# --- Cross-lineage shared ---------------------------------------------------

describe("Cross-lineage shared", {
  it("Case 1: Myeloid + Lymphoid (Mono, B)", {
    sig <- c("predicted.celltype.l1.Mono", "predicted.celltype.l1.B")
    expected <- "Cross-lineage shared"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 2: Myeloid + T-cell (DC, CD4_T)", {
    sig <- c("predicted.celltype.l1.DC", "predicted.celltype.l1.CD4_T")
    expected <- "Cross-lineage shared"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 13: All 8 cell types significant", {
    sig <- c(
      "predicted.celltype.l1.Mono",
      "predicted.celltype.l1.DC",
      "predicted.celltype.l1.NK",
      "predicted.celltype.l1.B",
      "predicted.celltype.l1.CD4_T",
      "predicted.celltype.l1.CD8_T",
      "predicted.celltype.l1.other_T",
      "predicted.celltype.l1.other"
    )
    expected <- "Cross-lineage shared"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })
})


# --- Likely shared but underpowered ------------------------------------------

describe("Likely shared but underpowered", {
  it("Case 3: PBMC only", {
    sig <- "predicted.celltype.l1.PBMC"
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 4: Mono only, gray zone in B (LFSR = 0.15)", {
    sig <- "predicted.celltype.l1.Mono"
    lfsr <- c("predicted.celltype.l1.B" = 0.15)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig, lfsr)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 8: Mono only (lineage-specific candidate), gray zone in CD4_T (LFSR = 0.3)", {
    sig <- "predicted.celltype.l1.Mono"
    lfsr <- c("predicted.celltype.l1.CD4_T" = 0.3)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig, lfsr)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 10: CD4_T only, gray zone in Mono (LFSR = 0.1)", {
    sig <- "predicted.celltype.l1.CD4_T"
    lfsr <- c("predicted.celltype.l1.Mono" = 0.1)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig, lfsr)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 14: NK only, gray zone in Mono (LFSR = 0.4)", {
    sig <- "predicted.celltype.l1.NK"
    lfsr <- c("predicted.celltype.l1.Mono" = 0.4)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig, lfsr)
    expect_equal(cpp_result, oracle_result)
  })
})


# --- Lineage-specific -------------------------------------------------------

describe("Lineage-specific", {
  it("Case 6: Mono + DC (myeloid only)", {
    sig <- c("predicted.celltype.l1.Mono", "predicted.celltype.l1.DC")
    expected <- "Lineage-specific"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 7: B + NK (lymphoid only)", {
    sig <- c("predicted.celltype.l1.B", "predicted.celltype.l1.NK")
    expected <- "Lineage-specific"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 18: other + Mono (other excluded from lineage sets)", {
    # sig_non_pbmc = {other, Mono}: 2 CTs, skip single-CT branch.
    # myeloid_sig = {Mono}, lymphoid_sig = {}
    # has_myeloid && !has_lymphoid -> Lineage-specific (no LFSR to trigger gray zone)
    sig <- c("predicted.celltype.l1.other", "predicted.celltype.l1.Mono")
    expected <- "Lineage-specific"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })
})


# --- T-cell-specific --------------------------------------------------------

describe("T-cell-specific", {
  it("Case 9: CD4_T + CD8_T", {
    sig <- c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T")
    expected <- "T-cell-specific"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 11: All T-cells (CD4_T + CD8_T + other_T)", {
    sig <- c(
      "predicted.celltype.l1.CD4_T",
      "predicted.celltype.l1.CD8_T",
      "predicted.celltype.l1.other_T"
    )
    expected <- "T-cell-specific"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })
})


# --- Single cell-type -------------------------------------------------------

describe("Single cell-type", {
  it("Case 5: Mono only, no gray zone (B LFSR = 0.8)", {
    sig <- "predicted.celltype.l1.Mono"
    lfsr <- c("predicted.celltype.l1.B" = 0.8)
    expected <- "Single cell-type"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig, lfsr)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 15: NK only, null in Mono (LFSR = 0.9)", {
    sig <- "predicted.celltype.l1.NK"
    lfsr <- c("predicted.celltype.l1.Mono" = 0.9)
    expected <- "Single cell-type"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig, lfsr)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 16: PBMC + Mono (PBMC stripped, leaving 1 CT)", {
    # PBMC is removed from sig_non_pbmc, leaving only Mono.
    # 1 CT, no LFSR provided -> Single cell-type.
    sig <- c("predicted.celltype.l1.PBMC", "predicted.celltype.l1.Mono")
    expected <- "Single cell-type"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })

  it("Case 17: other only (not in myeloid or lymphoid)", {
    # sig_non_pbmc = {other}: 1 CT. No LFSR -> Single cell-type.
    sig <- "predicted.celltype.l1.other"
    expected <- "Single cell-type"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })
})


# --- No significance --------------------------------------------------------

describe("No significance", {
  it("Case 12: No significant cell types", {
    sig <- character(0)
    expected <- "No significance"

    oracle_result <- oracle_categorize_specificity(sig)
    expect_equal(oracle_result, expected)

    cpp_result <- call_cpp(sig)
    expect_equal(cpp_result, oracle_result)
  })
})

# --- LFSR boundary tests ---------------------------------------------------

describe("LFSR boundary tests", {
  it("LFSR exactly 0.05 triggers gray zone (single CT branch)", {
    sig <- c("predicted.celltype.l1.Mono")
    lfsr <- c("predicted.celltype.l1.B" = 0.05)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })

  it("LFSR exactly 0.5 does NOT trigger gray zone", {
    sig <- c("predicted.celltype.l1.Mono")
    lfsr <- c("predicted.celltype.l1.B" = 0.5)
    expected <- "Single cell-type"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })

  it("LFSR just below 0.5 triggers gray zone", {
    sig <- c("predicted.celltype.l1.Mono")
    lfsr <- c("predicted.celltype.l1.B" = 0.4999)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })

  it("NaN LFSR does not trigger gray zone", {
    sig <- c("predicted.celltype.l1.Mono")
    lfsr <- c("predicted.celltype.l1.B" = NaN)
    expected <- "Single cell-type"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })
})

# --- T-cell underpowered branch --------------------------------------------

describe("T-cell underpowered branch", {
  it("T-cell with gray zone in myeloid → underpowered", {
    sig <- c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T")
    lfsr <- c("predicted.celltype.l1.Mono" = 0.15)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })

  it("T-cell with gray zone in non-T lymphoid → underpowered", {
    sig <- c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T")
    lfsr <- c("predicted.celltype.l1.B" = 0.2)
    expected <- "Likely shared but underpowered"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })

  it("T-cell with other in gray zone → stays T-cell-specific (other ignored)", {
    sig <- c("predicted.celltype.l1.CD4_T", "predicted.celltype.l1.CD8_T")
    lfsr <- c("predicted.celltype.l1.other" = 0.1)
    expected <- "T-cell-specific"

    oracle_result <- oracle_categorize_specificity(sig, lfsr_values = lfsr)
    expect_equal(oracle_result, expected)
    expect_equal(call_cpp(sig, lfsr), oracle_result)
  })
})
