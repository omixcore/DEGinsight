#' Example datasets bundled with DEGinsight
#'
#' @description
#' Simulated data mimicking a TNF-stimulation bulk RNA-seq experiment,
#' used in examples, tests and the vignette.
#'
#' \describe{
#'   \item{example_deg}{Simulated DEG table: `symbol`, `log2FoldChange`,
#'     `padj` for ~2,020 genes.}
#'   \item{example_gene_sets}{Small named-list gene-set collection
#'     (TNF/NF-kB, IL-17, apoptosis, housekeeping sets).}
#'   \item{example_ora}{Pre-computed ORA result on `example_deg`
#'     (via \code{\link{run_ora}}), so examples run without internet access.}
#' }
#'
#' @docType data
#' @name example_data
#' @aliases example_deg example_gene_sets example_ora
#' @usage data(example_data)
#' @format Three objects (see description).
#' @keywords datasets
NULL
