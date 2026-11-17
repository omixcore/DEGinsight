source("benchmark/common.R")
cfg <- model_config(Sys.getenv("ABLATION_MODEL", "deepseek-flash"))
if (!nzchar(cfg$api_key)) stop("Set credentials for the ablation model.")
datasets <- strsplit(Sys.getenv("ABLATION_DATASETS",
  "rapamycin_nk,imatinib_k562,dexamethasone_organoid"), ",")[[1]]
grid <- expand.grid(few_shot = c(FALSE, TRUE), top_n = c(5L, 10L, 20L, 30L),
                    tissue = c(FALSE, TRUE))
failures <- character()
for (id in datasets) {
  for (i in seq_len(nrow(grid))) {
    message(id, " / configuration ", i)
    tryCatch(run_cell(id, c(cfg, as.list(grid[i, ])), "ablation"),
      error = function(e) {
        failures <<- c(failures, paste(id, i, conditionMessage(e)))
        message(tail(failures, 1L))
      })
  }
}
if (length(failures)) stop(length(failures), " ablation cells failed. Successful cells were saved.")
