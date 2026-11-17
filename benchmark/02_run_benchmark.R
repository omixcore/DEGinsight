source("benchmark/common.R")
models <- strsplit(Sys.getenv("BENCH_MODELS", "deepseek-flash,deepseek-v4-pro"), ",")[[1]]
configs <- lapply(models, model_config)
if (any(!vapply(configs, function(x) nzchar(x$api_key), logical(1)))) {
  stop("Set the credentials for every requested model in benchmark/models.csv.")
}
failures <- character()
for (id in GT$dataset_id) {
  for (cfg in configs) {
    message(id, " / ", cfg$model)
    tryCatch(run_cell(id, c(cfg, list(tissue = TRUE, top_n = 20L, few_shot = FALSE)), "results"),
      error = function(e) {
        failures <<- c(failures, paste(id, cfg$model, conditionMessage(e)))
        message(tail(failures, 1L))
      })
  }
}
if (length(failures)) stop(length(failures), " benchmark cells failed. Successful cells were saved.")
