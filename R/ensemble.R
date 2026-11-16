#' Ensemble interpretation across multiple LLMs
#'
#' Implements the "elite ensemble" strategy: query several frontier models
#' independently, validate each output against controlled vocabularies, then
#' merge hypotheses by majority vote over supporting pathways. Hypotheses
#' whose supporting-pathway sets are cited by a strict minority of models
#' are downgraded to `speculative`.
#'
#' @inheritParams interpret_enrichment
#' @param few_shot Logical; passed through to
#'   \code{\link{interpret_enrichment}} (ablation factor).
#' @param models Character vector (length >= 2) of model IDs.
#' @param min_agreement Minimum fraction of models that must support a
#'   pathway cluster for it to enter the consensus (default 0.5).
#'
#' @return A `"deginsight_ensemble"` object: list with `consensus`
#'   (merged hypotheses), `per_model` (individual
#'   `"deginsight_interpretation"` objects) and `agreement` (pairwise
#'   pathway-level agreement matrix).
#'
#' @export
interpret_ensemble <- function(enrichment,
                               contrast,
                               species,
                               tissue = NULL,
                               direction = c("mixed", "up", "down"),
                               top_n = 20,
                               models = c("deepseek-chat", "kimi-k2", "gpt-4o-mini"),
                               vocabulary = NULL,
                               min_agreement = 0.5,
                               ...) {
  if (length(models) < 2) {
    stop("Ensemble needs >= 2 models; use interpret_enrichment() for one.")
  }
  direction <- match.arg(direction)

  per_model <- lapply(models, function(m) {
    interpret_enrichment(enrichment, contrast, species, tissue, direction,
                         top_n, model = m, vocabulary = vocabulary, ...)
  })
  names(per_model) <- models

  # pathway-level agreement: fraction of models citing each pathway
  all_cited <- unique(unlist(lapply(per_model, function(x) x$validation$cited)))
  support <- vapply(all_cited, function(pth) {
    mean(vapply(per_model, function(x) pth %in% x$validation$cited, logical(1)))
  }, numeric(1))

  consensus_pathways <- names(support)[support >= min_agreement]

  consensus <- list(
    pathways = consensus_pathways,
    support = support,
    hallucination_rates = vapply(per_model,
                                 function(x) x$validation$hallucination_rate,
                                 numeric(1))
  )

  out <- list(consensus = consensus, per_model = per_model,
              agreement = agreement_matrix(per_model))
  class(out) <- "deginsight_ensemble"
  out
}

#' Jaccard-similarity majority vote over sets of cited pathways
#'
#' @param sets A named list of character vectors (model -> cited pathways).
#' @param min_agreement Minimum fraction of sets an element must appear in.
#' @return Character vector of consensus elements.
#' @export
majority_vote <- function(sets, min_agreement = 0.5) {
  elems <- unique(unlist(sets))
  freq <- vapply(elems, function(e) {
    mean(vapply(sets, function(s) e %in% s, logical(1)))
  }, numeric(1))
  elems[freq >= min_agreement]
}

#' @noRd
agreement_matrix <- function(per_model) {
  models <- names(per_model)
  m <- matrix(NA_real_, length(models), length(models),
              dimnames = list(models, models))
  for (i in models) for (j in models) {
    a <- per_model[[i]]$validation$cited
    b <- per_model[[j]]$validation$cited
    m[i, j] <- if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else 1
  }
  m
}

#' @export
print.deginsight_ensemble <- function(x, ...) {
  cat("DEGinsight ensemble (", length(x$per_model), " models)\n", sep = "")
  cat("  consensus pathways:", length(x$consensus$pathways), "\n")
  cat("  hallucination rates:",
      paste(sprintf("%s=%.0f%%", names(x$consensus$hallucination_rates),
                    100 * x$consensus$hallucination_rates), collapse = ", "),
      "\n")
  invisible(x)
}
