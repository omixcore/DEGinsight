#' Validate LLM pathway claims against controlled vocabularies
#'
#' Hallucination control is central to DEGinsight. This function checks every
#' pathway name cited by the model against (a) the pathway names actually
#' present in the input enrichment table and (b) an optional larger reference
#' vocabulary (e.g. all GO Biological Process terms).
#'
#' Pathways cited by the model that occur in \emph{neither} set are flagged
#' as hallucinations. The hallucination rate is a headline metric in the
#' DEGinsight benchmark.
#'
#' @param interpretation Raw parsed LLM output (a list with `hypotheses`
#'   and optionally `flagged_pathways`), or a
#'   `"deginsight_interpretation"` object.
#' @param reference_terms Character vector of accepted term names
#'   (case-insensitive matching).
#' @param input_pathways Character vector of pathway names from the input
#'   enrichment table (these are always accepted).
#'
#' @return A list with:
#' \describe{
#'   \item{cited}{All distinct pathway names cited by the model.}
#'   \item{hallucinated}{Cited names found in neither reference set.}
#'   \item{hallucination_rate}{Proportion of cited terms that are
#'     hallucinated (0-1).}
#'   \item{per_hypothesis}{Data.frame of per-hypothesis validity.}
#' }
#'
#' @export
validate_terms <- function(interpretation, reference_terms, input_pathways) {
  hyps <- interpretation$hypotheses %||% list()

  cited <- as.character(unlist(lapply(hyps, function(h) h$supporting_pathways)))
  cited <- trimws(cited[!is.na(cited)])
  cited <- cited[nzchar(cited) & !duplicated(tolower(cited))]
  input_match <- match(tolower(cited), tolower(input_pathways))
  grounded <- unique(input_pathways[stats::na.omit(input_match)])

  ok <- tolower(trimws(unique(c(reference_terms, input_pathways))))
  is_valid <- tolower(cited) %in% ok
  hallucinated <- cited[!is_valid]

  per_hyp <- if (length(hyps)) {
    do.call(rbind, lapply(hyps, function(h) {
      sp <- trimws(as.character(unlist(h$supporting_pathways)))
      sp <- sp[!is.na(sp) & nzchar(sp)]
      data.frame(
        hypothesis = h$name %||% NA_character_,
        n_cited = length(sp),
        n_valid = sum(tolower(sp) %in% ok),
        stringsAsFactors = FALSE
      )
    }))
  } else data.frame()

  list(
    cited = cited,
    grounded = grounded,
    unsupported_pathways = cited[is.na(input_match)],
    hallucinated = hallucinated,
    hallucination_rate = if (length(cited)) length(hallucinated) / length(cited) else 0,
    per_hypothesis = per_hyp
  )
}

#' Load the GO Biological Process vocabulary as a reference
#'
#' Requires \pkg{GO.db} (Bioconductor).
#'
#' @param ontology One of `"BP"`, `"MF"`, `"CC"`.
#' @return Character vector of GO term names.
#' @export
go_vocabulary <- function(ontology = c("BP", "MF", "CC")) {
  ontology <- match.arg(ontology)
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    stop("Install GO.db first: BiocManager::install('GO.db')", call. = FALSE)
  }
  terms <- as.data.frame(GO.db::GOTERM)
  unique(terms$Term[terms$Ontology == ontology])
}
