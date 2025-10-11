test_that("TCGA_methylation_expression_correlation() works", {
    
    ## returns a tibble when values_only is set to TRUE
    res <- TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                   gene = "TDRD1",
                                                   values_only = TRUE)
    expect_true(inherits(res, "data.frame"))
    
    ## results are filtered according to the minimum number of probe to include
    res <- TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                   gene = "TDRD1",
                                                   min_probe_number = 6,
                                                   values_only = TRUE)
    expect_true(all(res$probe_number >= 6)) 
    
    ## returns an error if not enough probes to evaluate the mean methylation
    expect_error(TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                         gene = "MAGEA1",
                                                         min_probe_number = 3,
                                                         values_only = TRUE), 
                 "probes")
    
    ## returns a plot by default 
    ## res <- TCGA_methylation_expression_correlation(tumor = "LUAD", 
    ##                                               gene = "TDRD1")
    ## expect_s3_class(res, "ggplot")
    ## vdiffr::expect_doppelganger("TCGA_methylation_expression_correlation_TDRD1_LUAD", fig = res)
})
