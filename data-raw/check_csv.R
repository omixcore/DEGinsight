d <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
              check.names = FALSE)
cat("ncol:", ncol(d), " nrow:", nrow(d), "\n")
cat("colnames:", paste(names(d), collapse = " | "), "\n")
print(d[, c("dataset_id", "gse", "treat_regex")])
