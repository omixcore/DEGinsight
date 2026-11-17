#' Render a reproducible HTML interpretation report
#'
#' Renders the bundled R Markdown template with the enrichment table and
#' interpretation results.
#'
#' @param ora Enrichment data.frame (from \code{\link{run_ora}}).
#' @param interpretation A `"deginsight_interpretation"` or
#'   `"deginsight_ensemble"` object.
#' @param output_file Output HTML path.
#' @param title Report title.
#' @return The output path (invisibly).
#'
#' @export
render_report <- function(ora, interpretation,
                          output_file = "DEGinsight_report.html",
                          title = "DEGinsight mechanism interpretation") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Install rmarkdown first: install.packages('rmarkdown')", call. = FALSE)
  }
  template <- system.file("templates", "report.Rmd", package = "DEGinsight")
  if (!dir.exists(dirname(output_file))) stop("Report output directory does not exist.")
  output_file <- file.path(normalizePath(dirname(output_file)), basename(output_file))
  rmarkdown::render(template,
                    params = list(ora = ora, interpretation = interpretation,
                                  title = title),
                    output_file = output_file,
                    envir = new.env(parent = globalenv()))
  invisible(output_file)
}
