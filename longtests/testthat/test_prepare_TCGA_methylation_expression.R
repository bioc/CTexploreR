test_that("prepare_TCGA_methylation_expression_() works", {

  ## returns a tibble
  ## returns values for all tumors when tumor parameter is set to "all"
  ## Peritumoral samples are not included by default
  res <- CTexploreR:::prepare_TCGA_methylation_expression(tumor = "all",
                                                          gene = "TDRD1")
  expect_true(inherits(res, "data.frame"))
  expect_true(all(unique(res$type) %in% c("SKCM", "LUAD", "LUSC", "COAD",
                                          "ESCA", "BRCA", "HNSC")))
  expect_true(!"Peritumoral" %in% res$tissue)

  # ## returns an error if no valid tumor type is entered.
  # expect_error(prepare_TCGA_methylation_expression(
  #   tumor = "xxx",
  #   gene = "TDRD1"),
  #   "No valid")

  # ## returns an error if no valid gene name entered
  # expect_error(prepare_TCGA_methylation_expression(
  #   tumor = "LUAD",
  #   gene = "xxx"),
  #   "No valid")

  # ## returns an error if more than one gene entered
  # expect_error(prepare_TCGA_methylation_expression(
  #   tumor = "LUAD",
  #   gene = c("MAGEA1", "MAGEA3")),
  #   "No valid")

  ## returns values for the specified tumor type
  ## Peritumoral samples are returned when include_normal_tissues is TRUE
  res <- CTexploreR:::prepare_TCGA_methylation_expression(tumor = "LUAD",
                                                          gene = "TDRD1",
                                                          include_normal_tissues = TRUE)
  expect_true(all(unique(res$type) == "LUAD"))
  expect_true("Peritumoral" %in% res$tissue)
})
