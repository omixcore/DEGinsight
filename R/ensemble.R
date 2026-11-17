#' Ensemble interpretation across multiple LLMs
#'
#' Computes agreement over cited pathways supplied in the prompt.
#' Mechanisms enter consensus only when models use the same normalized claim
#' and share grounded pathway evidence. Original evidence grades are retained.
#'
#' @inheritParams interpret_enrichment
#' @param few_shot Logical; passed through to
#'   \code{\link{interpret_enrichment}} (ablation factor).
#' @param models Character vector (length >= 2) of model IDs.
#' @param min_agreement Minimum fraction of models that must support a
#'   pathway to enter the consensus. NULL requires a strict majority.
#'
#' @return A `"deginsight_ensemble"` object: list with `consensus`
#'   (grounded pathways, supporting models, and shared mechanism claims), `per_model` (individual
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
                               min_agreement = NULL,
                               ..., few_shot = FALSE) {
  if (length(models) < 2L || anyNA(models) || any(!nzchar(models)) || anyDuplicated(models)) {
    stop("Ensemble needs at least two distinct, nonempty model IDs.")
  }
  if (!is.null(min_agreement) && (length(min_agreement) != 1L ||
      !is.finite(min_agreement) || min_agreement <= 0 || min_agreement > 1)) {
    stop("min_agreement must be in (0, 1].")
  }
  direction <- match.arg(direction)

  per_model <- lapply(models, function(m) {
    interpret_enrichment(enrichment, contrast, species, tissue, direction,
                         top_n, model = m, vocabulary = vocabulary, few_shot = few_shot, ...)
  })
  names(per_model) <- models

  # pathway-level agreement: fraction of models citing each pathway
  all_cited <- unique(as.character(unlist(lapply(per_model, function(x) x$validation$grounded))))
  support <- vapply(all_cited, function(pth) {
    mean(vapply(per_model, function(x) pth %in% x$validation$grounded, logical(1)))
  }, numeric(1))

  keep <- if (is.null(min_agreement)) support > 0.5 else support >= min_agreement
  consensus_pathways <- all_cited[keep]
  supporting_models <- lapply(stats::setNames(consensus_pathways, consensus_pathways), function(pth) {
    models[vapply(per_model, function(x) pth %in% x$validation$grounded, logical(1))]
  })

  consensus <- list(
    pathways = consensus_pathways,
    support = support,
    supporting_models = supporting_models,
    mechanisms = mechanism_consensus(per_model, min_agreement),
    hallucination_rates = vapply(per_model,
                                 function(x) x$validation$hallucination_rate,
                                 numeric(1))
  )

  out <- list(consensus = consensus, per_model = per_model,
              agreement = agreement_matrix(per_model))
  class(out) <- "deginsight_ensemble"
  out
}

#' @noRd
mechanism_consensus <- function(per_model, min_agreement) {
  models <- names(per_model)
  evidence <- lapply(models, function(model) {
    result <- per_model[[model]]
    lapply(result$hypotheses, function(h) {
      keys <- tolower(trimws(as.character(unlist(h$supporting_pathways))))
      grounded <- result$validation$grounded
      supporting <- grounded[tolower(grounded) %in% keys]
      list(model = model, claim = h$mechanism, pathways = supporting,
           drivers = intersect(unlist(h$driver_genes), result$input_genes),
           grade = h$evidence_grade)
    })
  })
  evidence <- do.call(c, evidence)
  evidence <- Filter(function(x) length(x$pathways) > 0L && nzchar(trimws(x$claim)), evidence)
  keys <- vapply(evidence, function(x) normalize_omics_text(x$claim), character(1))
  groups <- split(evidence, keys)
  threshold <- function(n) {
    if (is.null(min_agreement)) n > length(models) / 2 else n / length(models) >= min_agreement
  }
  agreements <- lapply(groups, function(group) {
    pathways <- unique(unlist(lapply(group, function(x) x$pathways)))
    support <- lapply(stats::setNames(pathways, pathways), function(p) {
      unique(vapply(Filter(function(x) p %in% x$pathways, group),
                    function(x) x$model, character(1)))
    })
    support <- Filter(function(x) threshold(length(x)), support)
    if (!length(support)) return(NULL)
    list(mechanism = group[[1L]]$claim, pathways = names(support),
         supporting_models = support,
         evidence = Filter(function(x) any(x$pathways %in% names(support)), group))
  })
  unname(Filter(Negate(is.null), agreements))
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
    a <- per_model[[i]]$validation$grounded
    b <- per_model[[j]]$validation$grounded
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
