# Run from the repository root: Rscript benchmark/run_all.R [--live] [--gsea]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% c("--live", "--gsea")),
          requireNamespace("limma", quietly = TRUE),
          requireNamespace("edgeR", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))
run <- function(command, args) {
  status <- system2(command, args)
  if (status != 0L) stop(command, " ", paste(args, collapse = " "), " failed with status ", status)
}
rscript <- file.path(R.home("bin"), "Rscript")
run("bash", "benchmark/00_fetch_geo.sh")
run("python3", "benchmark/01_convert_inputs.py")
run(rscript, "benchmark/01b_prepare_real.R")
run(rscript, "benchmark/04b_baseline_ora.R")
if ("--gsea" %in% args) run(rscript, "benchmark/04c_baseline_gsea.R")
if ("--live" %in% args) {
  stages <- c("02_run_benchmark.R", "03_run_ablation.R", "04_score.R")
  status <- vapply(file.path("benchmark", stages), function(stage) system2(rscript, stage), integer(1))
  if (any(status != 0L)) stop("Benchmark stages failed: ", paste(names(status)[status != 0L], collapse = ", "))
}
