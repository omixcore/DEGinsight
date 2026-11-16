lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
options(timeout = 600, Ncpus = 4)
repos <- c(CRAN = "https://cloud.r-project.org")
cran_pkgs <- c("httr2", "jsonlite", "digest", "testthat", "knitr", "rmarkdown")
install.packages(cran_pkgs, repos = repos, lib = lib)
cat("CRAN done\n")
install.packages("BiocManager", repos = repos, lib = lib)
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE,
                     lib = lib)
cat("Bioconductor done\n")
for (p in c(cran_pkgs, "GEOquery", "limma")) {
  cat(p, requireNamespace(p, quietly = TRUE), "\n")
}
