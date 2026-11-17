#' Over-representation analysis (ORA) of DEGs
#'
#' A lightweight, dependency-free ORA implementation (hypergeometric test +
#' BH adjustment). Accepts any gene-set collection in GMT-like named-list
#' form, so results stay reproducible without external databases.
#'
#' For full-featured enrichment you can instead pass results from
#' \pkg{clusterProfiler} after mapping its columns to \code{\link{interpret_enrichment}} —
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
#'   adjusted p-value. Eligible zero-overlap pathways contribute to adjustment. Class `"deginsight_ora"`.
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
  size <- c(min_size, max_size)
  if (length(size) != 2L || any(!is.finite(size)) || any(size != floor(size)) ||
      min_size < 1 || max_size < min_size) stop("Invalid pathway size limits.")
  if (length(p_cutoff) != 1L || !is.finite(p_cutoff) || p_cutoff < 0 || p_cutoff > 1) {
    stop("p_cutoff must be in [0, 1].")
  }
  universe <- unique(as.character(universe))
  universe <- universe[!is.na(universe) & nzchar(universe)]
  if (!length(universe)) stop("The gene universe is empty.")
  query <- intersect(as.character(genes), universe)
  if (length(genes) && !length(query)) stop("No query genes overlap the universe.")
  genes <- query
  pathway_names <- names(gene_sets)
  if (!is.list(gene_sets) || (length(gene_sets) && is.null(pathway_names)) ||
      anyNA(pathway_names) || any(!nzchar(pathway_names)) || anyDuplicated(pathway_names)) {
    stop("gene_sets must be a list with unique, nonempty pathway names.")
  }
  N <- length(universe)
  n <- length(genes)

  res <- lapply(names(gene_sets), function(nm) {
    gs <- intersect(gene_sets[[nm]], universe)
    K <- length(gs)
    if (K < min_size || K > max_size) return(NULL)
    overlap <- intersect(genes, gs)
    k <- length(overlap)
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
    res <- data.frame(pathway = character(), overlap = integer(), size = integer(),
                      p_value = numeric(), genes = character())
  }
  res$p_adjust <- stats::p.adjust(res$p_value, method = p_adjust_method)
  res <- res[order(res$p_adjust, res$p_value, res$pathway, method = "radix"), c("pathway", "overlap", "size",
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

#' Preranked gene-set enrichment analysis
#'
#' Runs `fgsea::fgseaMultilevel()` on signed statistics for all tested genes.
#' Positive and negative normalized enrichment scores indicate concentration
#' at the top and bottom of the ranking, respectively.
#'
#' @param ranks Named numeric vector of finite, signed gene-level statistics.
#'   Use unique gene identifiers and avoid ties where possible.
#' @param gene_sets Named list of pathway gene vectors, as for
#'   \code{\link{run_ora}}.
#' @param min_size,max_size Pathway size filters within the ranked universe.
#' @return A data.frame with pathway, p_value, p_adjust, nes, size, and genes
#'   (comma-separated leading-edge genes).
#' @export
run_gsea <- function(ranks, gene_sets, min_size = 5L, max_size = 500L) {
  ids <- names(ranks)
  if (!is.numeric(ranks) || !length(ranks) || is.null(ids) ||
      any(!is.finite(ranks)) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("ranks must be a finite numeric vector with unique, nonempty gene names.")
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Install fgsea with BiocManager::install('fgsea').", call. = FALSE)
  }
  result <- fgsea::fgseaMultilevel(pathways = gene_sets, stats = ranks,
                                   minSize = min_size, maxSize = max_size, nproc = 1L)
  result <- result[is.finite(result$padj), ]
  out <- data.frame(pathway = result$pathway, p_value = result$pval,
                    p_adjust = result$padj, nes = result$NES, size = result$size,
                    genes = vapply(result$leadingEdge, paste, collapse = ",", FUN.VALUE = ""))
  out[order(out$p_adjust, out$p_value, out$pathway, method = "radix"), ]
}
