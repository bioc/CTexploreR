#' Check spelling of entered variables
#'
#' @noRd
#'
#' @description
#'
#' Checks the spelling of a vector of entered variable(s) comparing it
#' to a vector of valid names, and removes the ones that are absent
#' from the vector of valid names.
#'
#' @param variable `character()` containing the names of variables to
#'     check.
#'
#' @param valid_vector `character()` with valid variable names.
#'
#' @return A character with valid variables.
#'
#' @examples
#' CTexploreR:::check_names(
#'     variable = c("Ovarian", "leukemia", "wrong_name"),
#'     valid_vector = c("ovarian", "leukemia")
#' )
check_names <- function(variable, valid_vector) {
    in_valid <- variable %in% valid_vector
    if (!all(in_valid)) {
        msg <- paste0(
            sum(!in_valid), " out of ",
            length(in_valid), " names invalid: ",
            paste(variable[!in_valid], collapse = ", "),
            ".")
        warning(paste(strwrap(msg), collapse = "\n"),
            "\nSee the manual page for valid types.",
            call. = FALSE)
    }
    variable[in_valid]
}


#' Subset databases
#'
#' @noRd
#'
#' @description
#'
#' Check the presence of the genes in the database then subsets the database
#' to only keep these genes' data.
#'
#' @param variable `character()` containing the names genes to keep in the
#' data. The default value, `NULL`, takes CT (specific) genes
#'
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#' 
#' @param data `Summarized Experiment` or `SingleCellExperiment` object
#' with valid variable names.
#'
#' @return A `Summarized Experiment` or `SingleCellExperiment` object with
#' only the variables data
#'
#' @examples
#'
#'   CTexploreR:::subset_database(variable = "MAGEA1", data = CTdata::GTEX_data())
subset_database <- function(variable = NULL, data, include_CTP = FALSE) {
  CT_genes <- CTdata::CT_genes()
  if (is.null(variable) & include_CTP) 
    variable <- CT_genes$external_gene_name
  if (is.null(variable) & !include_CTP) 
    variable <- CT_genes[CT_genes$CT_gene_type == "CT_gene", 
                         ]$external_gene_name
  valid_gene_names <- unique(rowData(data)$external_gene_name)
  genes <- check_names(variable, valid_gene_names)
  database_subseted <- data[rowData(data)$external_gene_name %in% genes, ]
  return(database_subseted)
}

#' Prepare methylation and expression data of a gene in TCGA tumors
#'
#' @noRd
#'
#' @description
#'
#' Creates a Dataframe giving for each TCGA sample, the methylation level of
#' a gene (mean methylation of probes located in its promoter) and the
#' expression level of the gene (TPM value).
#'
#' @param gene `character` selected CT gene.
#'
#' @param tumor `character` defining the TCGA tumor type. Can be one
#'     of "SKCM", "LUAD", "LUSC", "COAD", "ESCA", "BRCA", "HNSC", or
#'     "all" (default).
#'
#' @param nt_up `numeric(1)` indicating the number of nucleotides
#'     upstream the TSS to define the promoter region (1000 by
#'     default)
#'
#' @param nt_down `numeric(1)` indicating the number of nucleotides
#'     downstream the TSS to define the promoter region (200 by
#'     default)
#'
#' @param include_normal_tissues `logical(1)`. If `TRUE`,
#'     the function will include normal peritumoral tissues in addition
#'     to tumoral samples. Default is `FALSE`.
#'
#' @return a Dataframe giving for each TCGA sample, the methylation level of
#' a gene (mean methylation of probes located in its promoter) and the
#' expression level of the gene (TPM value). The number of probes used to
#' estimate the methylation level is also reported.
#'
#' @importFrom BiocGenerics intersect
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#' @importFrom tibble tibble
#' @importFrom dplyr filter left_join select mutate
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   CTexploreR:::prepare_TCGA_methylation_expression("LUAD", gene = "TDRD1")
#' }
prepare_TCGA_methylation_expression <- function(tumor = "all",
                                                gene = NULL,
                                                nt_up = NULL,
                                                nt_down = NULL,
                                                include_normal_tissues = FALSE){
  suppressMessages({
    all_genes <- CTdata::all_genes()
    TPM <- CTdata::TCGA_TPM()
    met <- CTdata::TCGA_methylation()
  })

  TPM$type <- sub("TCGA-", "", TPM$project_id)
  met$type <- sub("TCGA-", "", met$project_id)
  valid_tumor <- c(unique(TPM$type), "all")
  type <- check_names(variable = tumor, valid_vector = valid_tumor)
  stopifnot("No valid tumor type entered" = length(type) > 0)
  if (!"all" %in% type) {
    TPM <- TPM[, TPM$type %in% type]
    met <- met[, met$type %in% type]
  }

  valid_gene_names <- all_genes$external_gene_name
  valid_gene_names <- valid_gene_names[valid_gene_names %in%
                                         rowData(TPM)$external_gene_name]
  gene <- check_names(gene, valid_gene_names)
  stopifnot("No valid gene name" = length(gene) == 1)
  TPM <- TPM[rowData(TPM)$external_gene_name %in% gene, ]

  ## Rm duplicated samples
  TPM <- TPM[, !duplicated(TPM$sample)]
  met <- met[, !duplicated(met$sample)]

  ## keep tumors for which expression and methylation data are available
  samples <- intersect(TPM$sample, met$sample)
  TPM <- TPM[, TPM$sample %in% samples]
  met <- met[, met$sample %in% samples]

  if (is.null(nt_up)) nt_up <- 1000
  if (is.null(nt_down)) nt_down <- 200

  ## Create a Grange corresponding to promoter region
  ## Calculates mean methylation value of promoter probe(s) in each sample
  gene_promoter_gr <- makeGRangesFromDataFrame(
    all_genes |>
      dplyr::filter(.data$external_gene_name == gene) |>
      dplyr::select(
        "ensembl_gene_id", "external_gene_name",
        "external_transcript_name", "chr", "strand",
        "transcription_start_site") |>
      mutate(chr = paste0("chr", .data$chr)) |>
      mutate(strand = ifelse(.data$strand == 1, "+", "-")) |>
      mutate(start = case_when(
        strand == "+" ~
          .data$transcription_start_site - nt_up,
        strand == "-" ~
          .data$transcription_start_site - nt_down)) |>
      mutate(stop = case_when(
        strand == "+" ~
          .data$transcription_start_site + nt_down,
        strand == "-" ~
          .data$transcription_start_site + nt_up)),
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "stop")

  met_roi <- subsetByOverlaps(met, gene_promoter_gr)
  met_mean <- colMeans(assay(met_roi), na.rm = TRUE)
  probe_number <- colSums(!is.na(assay(met_roi)))

  ## Keep prefix of TCGA sample names (TCGA-XX-XXXX-XXX) to join
  ## expression and methylation data
  names(met_mean) <- substr(names(met_mean), 1, 16)
  colnames(TPM) <- substr(colnames(TPM), 1, 16)

  methylation_expression <-
    suppressMessages(left_join(
      tibble(sample = names(met_mean),
             met = met_mean,
             probe_number = probe_number),
      tibble(sample = colnames(TPM),
             tissue = ifelse(TPM$shortLetterCode == "NT",
                             "Peritumoral", "Tumor"),
             type = TPM$type,
             TPM = as.vector(assay(TPM)))))

  if (!include_normal_tissues) methylation_expression <-
    methylation_expression[methylation_expression$tissue != "Peritumoral", ]


  methylation_expression <- dplyr::select(methylation_expression, "sample",
                                          "tissue", "type", "probe_number",
                                          "met", "TPM")

}






#' Determine font size
#'
#' @description
#'
#' Gives the fontsize to use for the heatmap based on the matrix's dimension.
#'
#' @param matrix `matrix` containing the data to visualise
#'
#' @return A logical number that is the fontsize to use
#'
#' @examples
#' CTexploreR:::set_fontsize(matrix(1:3, 9,8))
set_fontsize <- function(matrix) {
  if (dim(matrix)[1] > 140) return(3)
  if (dim(matrix)[1] > 100 & dim(matrix)[1] <= 140) return(4)
  if (dim(matrix)[1] > 50 & dim(matrix)[1] <= 100) return(5)
  if (dim(matrix)[1] > 20 & dim(matrix)[1] <= 50) return(6)
  if (dim(matrix)[1] <= 20) return(8)
}




#' CT genes description table
#'
#' @description
#'
#' Cancer-Testis (CT) genes description, imported from `CTdata`
#'
#' @format
#'
#' A `tibble` object with 280 rows and 47 columns.
#'
#' - Rows correspond to CT genes
#'
#' - Columns give CT genes characteristics
#'
#' @return A tibble of all 280 CT and CTP genes with their characteristics
#'
#' @details
#'
#' See `CTdata::CT_genes` documentation for details
#'
#' @source
#'
#' See `scripts/make_CT_genes.R` in `CTdata` for details on how this
#' list of curated CT genes was created.
#'
#' @export
#'
#' @name CT_genes
#'
#' @examples
#' CT_genes
#'
#' @docType data
CT_genes <- CTdata::CT_genes()


#' All genes description table
#'
#' @description
#'
#' All genes description, imported from `CTdata`
#'
#' @format
#'
#' A `tibble` object with 24488 rows and 47 columns.
#'
#' - Rows correspond to CT genes
#'
#' - Columns give CT genes characteristics
#'
#' @return A tibble of all 24 488 genes with their characteristics
#'
#' @details
#'
#' See `CTdata::all_genes` documentation for details
#'
#' @source
#'
#' See `scripts/make_all_genes.R` in `CTdata` for details on how this
#' list was created.
#'
#' @export
#'
#' @name all_genes
#'
#' @examples
#' all_genes
#'
#' @docType data
all_genes <- CTdata::all_genes()


## Heatmaps legend parameters

legends_param <- list(
  labels_gp = gpar(col = "black", fontsize = 6),
  title_gp = gpar(col = "black", fontsize = 6),
  row_names_gp = gpar(fontsize = 4),
  annotation_name_side = "left")

## Color vectors for heatmaps

legend_colors <- c(
    "#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
    "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F",
    "#9E0142")


CCLE_colors <- c(
    "lung" = "seagreen3", "skin" = "red3",
    "bile_duct" = "mediumpurple1", "bladder" = "mistyrose2",
    "colorectal" = "plum", "lymphoma" = "steelblue1",
    "uterine" = "darkorange4", "myeloma" = "turquoise3",
    "kidney" = "thistle4",
    "pancreatic" = "darkmagenta", "brain" = "palegreen2",
    "gastric" = "wheat3", "breast" = "midnightblue",
    "bone" = "sienna1", "head_and_neck" = "deeppink2",
    "ovarian" = "tan3", "sarcoma" = "lightcoral",
    "leukemia" = "steelblue4", "esophageal" = "khaki",
    "neuroblastoma" = "olivedrab1")


DAC_colors <- c(
    "B2-1" = "olivedrab2", "HCT116" = "lightcoral",
    "HEK293T" = "seagreen3", "HMLER" = "mediumpurple1",
    "IMR5-75" = "deeppink2", "NCH1681" = "steelblue2",
    "NCH612" = "red3", "TS603" = "darkmagenta")


TCGA_colors <- c(
    "BRCA" = "midnightblue", "COAD" = "darkorchid2",
    "ESCA" = "gold", "HNSC" = "deeppink2",
    "LUAD" = "seagreen", "LUSC" = "seagreen3",
    "SKCM" = "red3")

testis_colors <- c(
    "SSC" = "floralwhite", "Spermatogonia" = "moccasin",
    "Early_spermatocyte" = "gold",
    "Late_spermatocyte" = "orange",
    "Round_spermatid" = "red2",
    "Elongated_spermatid" = "darkred",
    "Sperm1" = "violet", "Sperm2" = "purple",
    "Sertoli" = "gray",
    "Leydig" = "cadetblue2", "Myoid" = "springgreen3",
    "Macrophage" = "gray10",
    "Endothelial" = "steelblue")

cell_type_colors <- c(
    "Germ_cells" = "mediumvioletred",
    "Trophoblast_cells" = "brown4",
    "Adipocyte_cells" = "lightgoldenrod1",
    "Blood_immune_cells" = "ivory2",
    "Endocrine_cells" = "lightslateblue",
    "Endothelial_cells" = "plum2",
    "Glandular_cells" = "deepskyblue",
    "Glial_cells" = "sienna1",
    "Mesenchymal_cells" = "peachpuff",
    "Pigment_cells" = "palegreen2", "Myoid" = "gray54",
    "Specialized_epithelial_cells" = "indianred1")

sex_colors <- c("F" = "deeppink", "M"= "steelblue")

Fetal_cell_colors = c("F_PGC" = "pink",
                      "F_GC" = "pink3",
                      "F_oogonia" = "palevioletred3",
                      "F_oocyte" = "mediumorchid4",
                      "M_PGC" = "lightblue1",
                      "M_GC" = "steelblue3",
                      "M_pre_spermatogonia" = "blue")

oocytes_colors <- c("Growing oocytes" = "mediumpurple1",
                    "Fully grown oocytes" ="lightgreen",
                    "Metaphase I" = "mediumvioletred",
                    "Metaphase II" = "mediumpurple4")

hESC_colors <- c("H7" = "salmon",
                 "H9" ="steelblue",
                 "H1" = "gray",
                 "HUES64" = "lightgoldenrod1")

Petropoulos_colors <- c("Morula E3" = "lightblue1",
                 "Blastocyst E4" ="lightblue3",
                 "Blastocyst E5" = "dodgerblue1",
                 "Blastocyst E6" = "dodgerblue4",
                 "Blastocyst E7" = "navy")

genotype_colors <- c("XX" = "pink", "XY"= "blue", "NA" = "gray", 
                     "X" = "pink4", "Y" = "blue4", "XXX"= "deeppink", "NA" = "gray")

fetal_time_colors <- colorRamp2(
  seq(6, 21, length = 6),                 
  c( "cyan", "cyan3", "dodgerblue", "dodgerblue3", "mediumpurple", "magenta4" ))

fetal_type_colors <- c("FGC" = "lightgreen", "Soma" = "gray")

embryo_stage_col <- c("GV Oocyte" = "mediumvioletred", "MII Oocyte" = "brown4", 
               "Sperm" = "palegreen2", "Zygote" = "lightgoldenrod1", 
               "2-cell" = "ivory2", "4-cell" = "lightslateblue", 
               "8-cell" = "plum2", "Morula" = "deepskyblue",
               "Blastocyst" = "sienna1",  "Post-implantation" = "peachpuff")

