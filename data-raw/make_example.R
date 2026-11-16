# Build the bundled example datasets (run from package root)
# Requires: source the package functions first.

source("R/cache.R")
source("R/enrichment.R")

set.seed(42)

# --- tiny built-in gene-set collection (TNF/NF-kB-flavoured) --------------
gene_pool <- paste0("GENE", sprintf("%04d", 1:2000))
example_gene_sets <- list(
  "TNF signaling via NF-kB" = c("TNF", "TNFR1", "NFKB1", "RELA", "IKBKB",
                                "CXCL8", "IL6", "ICAM1", "TRAF2", "BIRC3"),
  "IL-17 signaling" = c("IL17A", "IL17RA", "CXCL8", "CCL20", "MMP9",
                        "S100A8", "S100A9"),
  "Apoptosis" = c("CASP3", "CASP8", "BAX", "BCL2", "FAS", "TRADD"),
  "Ribosome biogenesis" = paste0("GENE", sprintf("%04d", 101:140)),
  "Oxidative phosphorylation" = paste0("GENE", sprintf("%04d", 201:260))
)
universe <- unique(c(unlist(example_gene_sets), gene_pool))

# --- simulated DEG table: TNF-stimulated vs control ------------------------
deg_genes <- c("TNF", "NFKB1", "RELA", "IKBKB", "CXCL8", "IL6", "ICAM1",
               "TRAF2", "BIRC3", "IL17A", "CCL20", "S100A8")
example_deg <- data.frame(
  symbol = universe,
  log2FoldChange = c(stats::rnorm(length(deg_genes), 2.2, 0.5),
                     stats::rnorm(length(universe) - length(deg_genes), 0, 0.6)),
  padj = c(stats::runif(length(deg_genes), 1e-8, 0.01),
           stats::runif(length(universe) - length(deg_genes), 0.05, 1)),
  stringsAsFactors = FALSE
)
example_deg$symbol[seq_along(deg_genes)] <- deg_genes

# --- pre-computed ORA result for offline examples/tests --------------------
sig <- example_deg$padj < 0.05 & abs(example_deg$log2FoldChange) > 1
example_ora <- run_ora(example_deg$symbol[sig], example_deg$symbol,
                       example_gene_sets, min_size = 3)

save(example_deg, example_gene_sets, example_ora,
     file = "data/example_data.rda", version = 2)
cat("Wrote data/example_data.rda\n")
