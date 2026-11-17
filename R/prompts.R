#' Build an ontology-constrained interpretation prompt
#'
#' Converts enrichment results into a structured prompt for mechanism
#' interpretation. The prompt design is a core methodological component of
#' DEGinsight:
#'
#' \itemize{
#'   \item \strong{Grounding}: the model is given the exact enriched term
#'     names, statistics, and overlapping genes — it is explicitly forbidden
#'     from inventing pathways.
#'   \item \strong{Context injection}: experimental contrast, species and
#'     tissue are provided to anchor biological reasoning.
#'   \item \strong{Evidence grading}: the model must grade each hypothesis
#'     (strong / moderate / speculative) based only on the supplied data.
#'   \item \strong{Structured output}: a JSON schema constrains the answer
#'     for downstream validation (\code{\link{validate_terms}}).
#' }
#'
#' @param enrichment A data.frame of enrichment results. Required columns:
#'   `pathway`, `p_adjust`; optional: `overlap`, `size`, `genes`.
#'   Output of \code{\link{run_ora}} or a `clusterProfiler` result coerced
#'   to the required columns.
#' @param contrast Character scalar describing the comparison, e.g.
#'   `"drug-treated vs control HepG2 cells"`.
#' @param species Character scalar, e.g. `"Homo sapiens"`.
#' @param tissue Optional tissue / cell-line context.
#' @param direction Character: are these up-regulated, down-regulated or
#'   mixed DEGs? One of `"up"`, `"down"`, `"mixed"`.
#' @param top_n Number of top pathways to include.
#' @param max_genes Maximum overlapping genes listed per pathway.
#' @param few_shot Logical; prepend one worked example (input -> expected
#'   JSON output) to the prompt. One of the ablation factors in the
#'   DEGinsight benchmark.
#'
#' @return A list with elements `system`, `prompt` and `schema`
#'   (the JSON response schema).
#'
#' @examples
#' data(example_data, package = "DEGinsight")
#' p <- build_interpretation_prompt(example_ora,
#'                                  contrast = "TNF-stimulated vs control",
#'                                  species = "Homo sapiens")
#' cat(substr(p$prompt, 1, 200))
#'
#' @export
build_interpretation_prompt <- function(enrichment,
                                        contrast,
                                        species,
                                        tissue = NULL,
                                        direction = c("mixed", "up", "down"),
                                        top_n = 20,
                                        max_genes = 30,
                                        few_shot = FALSE) {
  direction <- match.arg(direction)
  if (!is.data.frame(enrichment) || !all(c("pathway", "p_adjust") %in% names(enrichment)) ||
      !is.character(enrichment$pathway) || anyNA(enrichment$pathway) ||
      any(!nzchar(enrichment$pathway)) || !is.numeric(enrichment$p_adjust) ||
      any(!is.finite(enrichment$p_adjust) | enrichment$p_adjust < 0 | enrichment$p_adjust > 1)) {
    stop("Supply enrichment with pathway names and finite p_adjust values in [0, 1].")
  }
  if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1 || top_n != floor(top_n) ||
      length(max_genes) != 1L || !is.finite(max_genes) ||
      max_genes < 0 || max_genes != floor(max_genes)) stop("Invalid top_n or max_genes.")
  df <- enrichment
  if ("overlap" %in% names(df)) df <- df[is.na(df$overlap) | df$overlap > 0, , drop = FALSE]
  df <- utils::head(df[order(df$p_adjust), , drop = FALSE], top_n)

  shown_genes <- lapply(seq_len(nrow(df)), function(i) {
    if ("genes" %in% names(df) && !is.na(df$genes[i])) {
      utils::head(trimws(strsplit(as.character(df$genes[i]), ",")[[1]]), max_genes)
    } else character()
  })
  pathway_lines <- vapply(seq_len(nrow(df)), function(i) {
    genes <- paste(shown_genes[[i]], collapse = ", ")
    details <- if ("nes" %in% names(df)) {
      sprintf("NES = %+.2f, size = %s", df$nes[i], df$size[i])
    } else if (all(c("overlap", "size") %in% names(df))) {
      sprintf("overlap = %s/%s", df$overlap[i], df$size[i])
    } else "size and overlap unavailable"
    sprintf("%d. %s (adj.P = %.2e, %s)%s",
            i, df$pathway[i], df$p_adjust[i],
            details,
            if (nzchar(genes)) paste0("\n   Genes: ", genes) else "")
  }, character(1))

  system <- paste(
    "You are a careful molecular biologist specializing in transcriptomics.",
    "You interpret pathway enrichment results strictly based on the data",
    "provided. You NEVER invent pathway names, genes, or citations.",
    "When evidence is weak you say so explicitly."
  )

  prompt <- paste0(
    "## Experiment\n",
    "- Contrast: ", contrast, "\n",
    "- Species: ", species, "\n",
    if (!is.null(tissue)) paste0("- Tissue/cells: ", tissue, "\n") else "",
    "- Direction of DEGs analysed: ", direction, "-regulated\n\n",
    "## Enriched pathways (top ", nrow(df), " by adjusted p-value)\n",
    paste(pathway_lines, collapse = "\n"), "\n\n",
    if (few_shot) paste0(few_shot_example(), "\n\n") else "",
    "## Task\n",
    "1. Synthesize these pathways into 1-5 coherent mechanistic hypotheses ",
    "explaining the biological response in this contrast.\n",
    "2. For each hypothesis: name it, describe the mechanism, list the ",
    "supporting pathway names EXACTLY as written above, list key driver ",
    "genes (only from the genes shown), and grade the evidence as ",
    "'strong', 'moderate' or 'speculative'.\n",
    "3. Flag any pathway that looks like a false positive or artifact.\n",
    "4. Suggest one concrete follow-up experiment per hypothesis.\n",
    "NEVER invent pathway names, genes, or citations. ",
    "Use only evidence from this experiment. Enrichment and NES do not establish pathway activation or causality. ",
    "If no pathways were supplied, return empty hypotheses and flagged_pathways arrays."
  )

  schema <- list(
    hypotheses = paste(
      "array of objects: {name: string, mechanism: string,",
      "supporting_pathways: array of string, driver_genes: array of string,",
      "evidence_grade: 'strong'|'moderate'|'speculative',",
      "follow_up_experiment: string}"
    ),
    flagged_pathways = "array of objects: {pathway: string, reason: string}",
    overall_summary = "string"
  )

  list(system = system, prompt = prompt, schema = schema,
       input_pathways = df$pathway, input_genes = unique(unlist(shown_genes)))
}

#' Worked example used for few-shot prompting
#' @noRd
few_shot_example <- function() {
  paste0(
    "## Example (illustrating the expected reasoning and format)\n",
    "Input pathways:\n",
    "1. Example pathway A (adj.P = 1.2e-08, overlap = 18/120)\n   Genes: GENE_A, GENE_B\n",
    "2. Example pathway B (adj.P = 3.4e-06, overlap = 15/90)\n   Genes: GENE_C, GENE_D\n",
    "Expected output:\n",
    '```json',
    "\n{\n",
    '  "hypotheses": [{\n',
    '    "name": "A shared transcriptional response",\n',
    '    "mechanism": "Both gene sets are enriched and may reflect a shared response.",\n',
    '    "supporting_pathways": ["Example pathway A", "Example pathway B"],\n',
    '    "driver_genes": ["GENE_A", "GENE_C"],\n',
    '    "evidence_grade": "speculative",\n',
    '    "follow_up_experiment": "Measure GENE_A and GENE_C in independent samples."\n',
    "  }],\n",
    '  "flagged_pathways": [],\n',
    '  "overall_summary": "Two gene sets support a hypothesis of a shared response."\n',
    "}\n```"
  )
}
