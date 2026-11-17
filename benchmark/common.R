# Source the checkout under test, even if another DEGinsight is installed.
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
GT <- read.csv("benchmark/ground_truth.csv", fileEncoding = "UTF-8-BOM")
gene_sets <- read_gmt("benchmark/gene_sets.gmt")
MODELS <- read.csv("benchmark/models.csv", fileEncoding = "UTF-8-BOM")
BENCH <- Sys.getenv("DEGINSIGHT_BENCH_DIR", "benchmark/revised")
paths <- c("benchmark/ground_truth.csv", "benchmark/gene_sets.gmt", "benchmark/models.csv",
           "benchmark/sources.tsv", list.files("R", full.names = TRUE),
           list.files("benchmark", pattern = "\\.(R|py)$", full.names = TRUE))
protocol <- digest::digest(unname(tools::md5sum(paths)), algo = "sha256")
get_dataset <- function(id) {
  obj <- readRDS(file.path("benchmark/deg_tables", paste0(id, ".rds")))
  if (!identical(obj$pipeline, "limma-v2") ||
      !identical(obj$preparation_sha256, digest::digest(file = "benchmark/01b_prepare_real.R", algo = "sha256")) ||
      !identical(obj$sources_sha256, digest::digest(file = "benchmark/sources.tsv", algo = "sha256"))) {
    stop("Rebuild DEG tables with the current benchmark pipeline: ", id)
  }
  deg <- obj$deg
  ok <- !is.na(deg$symbol) & nzchar(deg$symbol) & is.finite(deg$padj) & is.finite(deg$log2FoldChange)
  deg <- deg[ok, ]
  # Select a representative independently of the DE statistic.
  deg <- deg[order(-deg$mean_expression, deg$gene_id, method = "radix"), ]
  deg <- deg[!duplicated(deg$symbol), ]
  obj$deg <- deg
  obj
}
benchmark_ora <- function(deg, direction = "mixed") {
  sel <- deg$padj < 0.05 & abs(deg$log2FoldChange) > 1
  if (direction == "up") sel <- sel & deg$log2FoldChange > 0
  if (direction == "down") sel <- sel & deg$log2FoldChange < 0
  genes <- deg$symbol[sel]
  ora <- run_ora(genes, deg$symbol, gene_sets)
  list(ora = ora[ora$overlap > 0, , drop = FALSE], n_degs = length(genes))
}
model_config <- function(id) {
  m <- MODELS[match(id, MODELS$model_id), ]
  if (nrow(m) != 1 || is.na(m$model_id)) stop("Unknown model: ", id)
  list(model = id, base_url = Sys.getenv(m$base_url_env, m$default_base_url),
       api_key = Sys.getenv(m$api_key_env))
}
save_run <- function(obj, kind, signature) {
  folder <- file.path(BENCH, kind)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(folder, paste0(signature, ".rds"))
  tmp <- tempfile(tmpdir = folder)
  on.exit(unlink(tmp))
  saveRDS(obj, tmp)
  if (!file.rename(tmp, out)) stop("Could not write ", out)
}
run_cell <- function(id, cfg, kind) {
  d <- GT[match(id, GT$dataset_id), ]
  if (is.na(d$dataset_id)) stop("Unknown dataset: ", id)
  ds <- get_dataset(id)
  b <- benchmark_ora(ds$deg)
  signature <- digest::digest(list(protocol, id, ds$deg, cfg[names(cfg) != "api_key"]),
                              algo = "sha256")
  path <- file.path(BENCH, kind, paste0(signature, ".rds"))
  if (file.exists(path)) return(invisible(NULL))
  started <- Sys.time()
  interp <- interpret_enrichment(b$ora, contrast = "treated vs control",
    species = "Homo sapiens", tissue = if (cfg$tissue) d$cell_type else NULL,
    model = cfg$model, base_url = cfg$base_url, api_key = cfg$api_key,
    top_n = cfg$top_n, few_shot = cfg$few_shot)
  save_run(list(dataset = id, model = cfg$model, protocol = protocol,
                config = cfg[names(cfg) != "api_key"], ora = b$ora, n_degs = b$n_degs,
                status = if (nrow(b$ora)) "completed" else "no_overlap",
                elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
                interpretation = interp), kind, signature)
}
read_runs <- function(kind) {
  rows <- lapply(list.files(file.path(BENCH, kind), "\\.rds$", full.names = TRUE), readRDS)
  rows <- Filter(function(x) identical(x$protocol, protocol), rows)
  rows <- Filter(function(x) {
    ds <- get_dataset(x$dataset)
    signature <- digest::digest(list(protocol, x$dataset, ds$deg, x$config), algo = "sha256")
    file.exists(file.path(BENCH, kind, paste0(signature, ".rds")))
  }, rows)
  rows
}
