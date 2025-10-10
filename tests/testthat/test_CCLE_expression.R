test_that("CCLE_expression() works", {
  
  ## returns a warning when an invalid gene is entered
  ## returns a matrix of double when values_only is set to TRUE
  ## works with only one gene in input
  ## tumor type is case insensitive
  ## selects the expected cell lines
  ## returns the expected matrix
  expect_warning(res <- CCLE_expression(c("MAGEA1", ""), 
                                        type = "sKIN", 
                                        values_only = TRUE), "names invalid")
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 1) 
  x <- colData(CCLE_data())
  exp_cells <- rownames(x[x$type == "Skin" , ])
  expect_equal(sort(colnames(res)), sort(exp_cells))
  #saveRDS(res, test_path("fixtures", "CCLE_expression_on_MAGE.rds"))
  exp_res <- readRDS(test_path("fixtures", "CCLE_expression_on_MAGE.rds"))
  expect_equal(exp_res, res)
  
  ## Test the "log_TPM" units argument
  res_in_log <- CCLE_expression(c("MAGEA1"), type = "sKIN", 
                                units = "log_TPM", 
                                values_only = TRUE)
  expect_equal(log1p(res), res_in_log) 
  
  ## returns a heatmap by default
  ## returns a warning when an invalid tumor type is entered
  ## expect_warning(res <- CCLE_expression(c("MAGEA1", "MAGEA3", "MAGEA4"), 
  ##                                         type = c("lung", "xxx")), 
  ##                  "names invalid")
  ##   expect_s4_class(res, "Heatmap")
  ##   vdiffr::expect_doppelganger("CCLE_expression_on_MAGE", fig = res)
  
  ## No valid tumor type returns an error
  expect_error(CCLE_expression(genes = "MAGEA1"), "No valid")
})



