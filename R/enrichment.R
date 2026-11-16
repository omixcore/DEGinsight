#' Over-representation analysis (ORA) of DEGs
#'
#' A lightweight, dependency-free ORA implementation (hypergeometric test +
#' BH adjustment). Accepts any gene-set collection in GMT-like named-list
#' form, so results stay reproducible without external databases.
#'
#' For full-featured enrichment you can instead pass results from
#' \pkg{clusterProfiler} directly to \code{\link{interpret_enrichment}} —
#' see the vignette.
#'
#' @param genes Character vector of significant genes (e.g. adjusted
#'   p < 0.05).
#' @param universe Character vector of the background gene set (all tested
#'   genes).
#' @param gene_sets Named list of character vectors: pathway name -> genes.
#'   See \code{\link{read_gmt}} for loading GMT files.
#' @param min_size,max_size Pathway size filters within the universe.
#' @param p_adjust_method Passed to \code{\link[stats]{p.adjust}}.
#' @param p_cutoff Keep pathways with adjusted p below this value (default
#'   keeps everything; filtering is usually done downstream).
#'
#' @return A data.frame with columns `pathway`, `overlap`, `size`,
#'   `p_value`, `p_adjust`, `genes` (comma-separated overlaps), ordered by
#'   adjusted p-value. Class `"deginsight_ora"`.
#'
#' @importFrom stats phyper p.adjust
#' @export
run_ora <- function(genes,
                    universe,
                    gene_sets,
                    min_size = 5,
                    max_size = 500,
                    p_adjust_method = "BH",
                    p_cutoff = 1) {
  genes <- intersect(unique(genes), universe)
  if (length(genes) == 0) stop("No query genes overlap the universe.")
  N <- length(universe)
  n <- length(genes)

  res <- lapply(names(gene_sets), function(nm) {
    gs <- intersect(gene_sets[[nm]], universe)
    K <- length(gs)
    if (K < min_size || K > max_size) return(NULL)
    overlap <- intersect(genes, gs)
    k <- length(overlap)
    if (k == 0) return(NULL)
    data.frame(
      pathway = nm,
      overlap = k,
      size = K,
      p_value = stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE),
      genes = paste(overlap, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, res)
  if (is.null(res) || nrow(res) == 0) {
    stop("No pathways passed the size filters / had any overlapping genes.")
  }
  res$p_adjust <- stats::p.adjust(res$p_value, method = p_adjust_method)
  res <- res[order(res$p_adjust), c("pathway", "overlap", "size",
                                    "p_value", "p_adjust", "genes")]
  rownames(res) <- NULL
  res <- res[res$p_adjust <= p_cutoff, ]
  class(res) <- c("deginsight_ora", "data.frame")
  res
}

#' Read a GMT file into a named list of gene sets
#'
#' @param path Path to a `.gmt` file.
#' @return Named list: gene set name -> character vector of genes.
#' @export
read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- lapply(strsplit(lines, "\t"), function(x) {
    if (length(x) < 3) return(NULL)
    x[-c(1, 2)]
  })
  names(out) <- vapply(strsplit(lines, "\t"), `[`, character(1), 1)
  out[!vapply(out, is.null, logical(1))]
}
