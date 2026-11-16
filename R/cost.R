#' Estimate the API cost of an interpretation run
#'
#' Computes a rough token count for the prompt and multiplies by per-model
#' prices, so users can budget before running large benchmarks.
#'
#' @param enrichment Enrichment results (used to build the prompt).
#' @param contrast,species,... Passed to
#'   \code{\link{build_interpretation_prompt}}.
#' @param models Character vector of model IDs to price.
#' @param prices Named list of per-model prices in USD per 1M tokens:
#'   each element `list(input = ..., output = ...)`. A small built-in table
#'   (\code{default_prices}) is provided; verify against current pricing.
#' @param expected_output_tokens Assumed completion length per call.
#'
#' @return A data.frame with per-model estimated cost in USD.
#'
#' @examples
#' data(example_ora, package = "DEGinsight")
#' estimate_cost(example_ora, contrast = "TNF vs control",
#'               species = "Homo sapiens",
#'               models = c("deepseek-chat", "gpt-4o-mini"))
#'
#' @export
estimate_cost <- function(enrichment,
                          contrast,
                          species,
                          ...,
                          models = c("deepseek-chat", "gpt-4o-mini"),
                          prices = default_prices(),
                          expected_output_tokens = 1500) {
  p <- build_interpretation_prompt(enrichment, contrast, species, ...)
  in_tokens <- nchar(p$prompt) %/% 4 + nchar(p$system) %/% 4  # ~4 chars/token

  do.call(rbind, lapply(models, function(m) {
    pr <- prices[[m]]
    if (is.null(pr)) {
      warning("No price for model '", m, "'; returning NA.", call. = FALSE)
      pr <- list(input = NA_real_, output = NA_real_)
    }
    data.frame(
      model = m,
      input_tokens = in_tokens,
      output_tokens = expected_output_tokens,
      cost_usd = in_tokens / 1e6 * pr$input +
        expected_output_tokens / 1e6 * pr$output,
      stringsAsFactors = FALSE
    )
  }))
}

#' @noRd
default_prices <- function() {
  # USD per 1M tokens, verified against provider pricing pages 2026-11-16.
  # DeepSeek publishes peak and off-peak rates (off-peak = 50% of peak);
  # the peak values below are the conservative upper bound. Peak hours are
  # 01:00-04:00 and 06:00-10:00 UTC on weekdays.
  list(
    "deepseek-flash"   = list(input = 0.30, output = 1.20),
    "deepseek-v4-pro"  = list(input = 1.32, output = 3.96),
    "deepseek-chat"    = list(input = 0.27, output = 1.10),
    "gpt-4o-mini"      = list(input = 0.15, output = 0.60),
    "gpt-4o"           = list(input = 2.50, output = 10.00),
    "kimi-k2"          = list(input = 0.60, output = 2.50),
    "claude-sonnet"    = list(input = 3.00, output = 15.00)
  )
}
