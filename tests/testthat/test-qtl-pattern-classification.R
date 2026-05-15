# Test all 25 QTL patterns via oracle functions
# Data-driven: define specs, loop over them

pattern_specs <- list(
  # 1. Local Cascade
  list(
    pattern = 1, mechanism = "Local Cascade",
    overlap = TRUE, caqtl_overlap = TRUE, caqtl_nonoverlap = FALSE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = TRUE, eqtl_linked_nonoverlap = FALSE
  ),
  # 2. Positional Cascade (overlap link, non-overlap caQTL)
  list(
    pattern = 2, mechanism = "Positional Cascade",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = TRUE, eqtl_linked_nonoverlap = FALSE
  ),
  # 3. Positional Cascade (overlap link)
  list(
    pattern = 3, mechanism = "Positional Cascade",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = TRUE, eqtl_linked_nonoverlap = FALSE
  ),
  # 4. Positional Cascade (non-overlap link)
  list(
    pattern = 4, mechanism = "Positional Cascade",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = TRUE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = TRUE
  ),
  # 5. Distal Cascade
  list(
    pattern = 5, mechanism = "Distal Cascade",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = TRUE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = TRUE
  ),
  # 6. caQTL + eQTL (overlap link, discordant)
  list(
    pattern = 6, mechanism = "caQTL + eQTL (No Link)",
    overlap = TRUE, caqtl_overlap = TRUE, caqtl_nonoverlap = FALSE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 7. caQTL + eQTL (overlap link, non-overlap caQTL, discordant)
  list(
    pattern = 7, mechanism = "caQTL + eQTL (No Link)",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 8. caQTL + eQTL (non-overlap link, discordant)
  list(
    pattern = 8, mechanism = "caQTL + eQTL (No Link)",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = TRUE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 9. caQTL + eQTL (no overlap, discordant)
  list(
    pattern = 9, mechanism = "caQTL + eQTL (No Link)",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = TRUE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 10. caQTL + eQTL (no link)
  list(
    pattern = 10, mechanism = "caQTL + eQTL (No Link)",
    overlap = TRUE, caqtl_overlap = TRUE, caqtl_nonoverlap = FALSE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 11. caQTL + eQTL (non-overlap caQTL, no link)
  list(
    pattern = 11, mechanism = "caQTL + eQTL (No Link)",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 12. caQTL + eQTL (no overlap, no link)
  list(
    pattern = 12, mechanism = "caQTL + eQTL (No Link)",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 13. Only caQTL (overlap link)
  list(
    pattern = 13, mechanism = "Only caQTL (With Link)",
    overlap = TRUE, caqtl_overlap = TRUE, caqtl_nonoverlap = FALSE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 14. Only caQTL (overlap link, non-overlap caQTL)
  list(
    pattern = 14, mechanism = "Only caQTL (With Link)",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 15. Only caQTL (non-overlap link)
  list(
    pattern = 15, mechanism = "Only caQTL (With Link)",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = TRUE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 16. Only caQTL (no overlap, link)
  list(
    pattern = 16, mechanism = "Only caQTL (With Link)",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = TRUE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 17. Only caQTL (overlap, no link)
  list(
    pattern = 17, mechanism = "Only caQTL (No Link)",
    overlap = TRUE, caqtl_overlap = TRUE, caqtl_nonoverlap = FALSE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 18. Only caQTL (non-overlap caQTL, no link)
  list(
    pattern = 18, mechanism = "Only caQTL (No Link)",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 19. Only caQTL (no overlap, no link)
  list(
    pattern = 19, mechanism = "Only caQTL (No Link)",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = TRUE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 20. Only eQTL (overlap link)
  list(
    pattern = 20, mechanism = "Only eQTL",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 21. Only eQTL (no link)
  list(
    pattern = 21, mechanism = "Only eQTL",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 22. Only eQTL (no overlap)
  list(
    pattern = 22, mechanism = "Only eQTL",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = TRUE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 23. No molQTL (overlap link)
  list(
    pattern = 23, mechanism = "No molQTL",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = TRUE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 24. No molQTL (no link)
  list(
    pattern = 24, mechanism = "No molQTL",
    overlap = TRUE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  ),
  # 25. No molQTL (no overlap)
  list(
    pattern = 25, mechanism = "No molQTL",
    overlap = FALSE, caqtl_overlap = FALSE, caqtl_nonoverlap = FALSE,
    link_overlap = FALSE, link_nonoverlap = FALSE,
    eqtl = FALSE, eqtl_linked_overlap = FALSE, eqtl_linked_nonoverlap = FALSE
  )
)

for (spec in pattern_specs) {
  test_that(paste0("Pattern ", spec$pattern, " - ", spec$mechanism), {
    result <- oracle_derive_pattern(
      has_overlap = spec$overlap,
      has_caqtl_overlap = spec$caqtl_overlap,
      has_caqtl_nonoverlap = spec$caqtl_nonoverlap,
      has_link_overlap = spec$link_overlap,
      has_link_nonoverlap = spec$link_nonoverlap,
      has_eqtl = spec$eqtl,
      eqtl_for_linked_overlap = spec$eqtl_linked_overlap,
      eqtl_for_linked_nonoverlap = spec$eqtl_linked_nonoverlap
    )
    expect_equal(result, spec$pattern)
    expect_equal(oracle_pattern_to_mechanism(result), spec$mechanism)
  })
}
