#' On-disk cache for LLM responses
#'
#' Creates (if needed) and returns the cache directory used by
#' \code{\link{call_llm}}. Location priority:
#' \enumerate{
#'   \item \code{dir} argument
#'   \item environment variable \code{DEGINSIGHT_CACHE_DIR}
#'   \item \code{tools::R_user_dir("DEGinsight", "cache")}
#' }
#'
#' @param dir Optional explicit directory.
#' @return The cache directory path (invisibly).
#' @export
DEGinsight_cache <- function(dir = NULL) {
  d <- dir %||% Sys.getenv("DEGINSIGHT_CACHE_DIR",
                           tools::R_user_dir("DEGinsight", "cache"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(d)
}

`%||%` <- function(x, y) {
  empty <- is.null(x) || length(x) == 0 ||
    (is.character(x) && length(x) == 1 && !nzchar(x))
  if (empty) y else x
}

#' @noRd
cache_path <- function(key) {
  file.path(DEGinsight_cache(), paste0(key, ".rds"))
}

#' @noRd
cache_get <- function(key) {
  p <- cache_path(key)
  if (file.exists(p)) readRDS(p) else NULL
}

#' @noRd
cache_set <- function(key, value) {
  saveRDS(value, cache_path(key))
  invisible(value)
}
