#' On-disk cache for LLM responses
#'
#' Creates (if needed) and returns the cache directory used by
#' \code{\link{call_llm}}. Location priority:
#' \enumerate{
#'   \item \code{dir} argument
#'   \item session option \code{DEGinsight.cache_dir}
#'   \item environment variable \code{DEGINSIGHT_CACHE_DIR}
#'   \item \code{tools::R_user_dir("DEGinsight", "cache")}
#' }
#'
#' @param dir Optional directory to use for the current R session.
#' @return The cache directory path (invisibly).
#' @export
DEGinsight_cache <- function(dir = NULL) {
  d <- dir %||% getOption("DEGinsight.cache_dir") %||%
    Sys.getenv("DEGINSIGHT_CACHE_DIR", tools::R_user_dir("DEGinsight", "cache"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Cannot create the DEGinsight cache directory.")
  d <- normalizePath(d)
  if (!is.null(dir)) options(DEGinsight.cache_dir = d)
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
  if (file.exists(p)) tryCatch(readRDS(p), error = function(e) NULL) else NULL
}

#' @noRd
cache_set <- function(key, value) {
  path <- cache_path(key)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(value, tmp)
  if (!file.exists(path) || is.null(cache_get(key))) {
    if (!suppressWarnings(file.rename(tmp, path)) && !file.exists(path)) {
      warning("Could not save the DEGinsight cache entry.")
    }
  }
  invisible(value)
}
