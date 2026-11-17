# DEGinsight

DEGinsight is an R package for interpreting bulk RNA-seq differential
expression with gene-set enrichment and large language models. It supplies
the model with pathway statistics and supporting genes, audits cited terms
and driver genes, and reports hypotheses and model agreement.

## Installation

```r
remotes::install_github("omixcore/DEGinsight")
```

R 4.1 or later is required. Install `fgsea` via BiocManager for preranked GSEA
and `rmarkdown` plus Pandoc for HTML reports.

## Example

```r
library(DEGinsight)
data(example_data)

# ORA runs locally. Use the genes tested in the DE model as the universe.
significant <- with(example_deg, padj < 0.05 & abs(log2FoldChange) > 1)
ora <- run_ora(example_deg$symbol[significant], example_deg$symbol,
               example_gene_sets)

# Set LLM_API_KEY and LLM_BASE_URL in the environment before calling a model.
result <- interpret_deg(
  example_deg, gene_col = "symbol", padj_col = "padj",
  logfc_col = "log2FoldChange", gene_sets = example_gene_sets,
  contrast = "treated vs control", species = "Homo sapiens",
  direction = "split"
)
render_report(result$up$ora, result$up$interpretation, "upregulated.html")
```

- `run_ora()` uses the hypergeometric test and adjusts over all size-eligible
  pathways. `interpret_deg(direction = "split")` runs separate up and down
  queries with a common tested-gene background.
- `run_gsea()` takes signed gene-level statistics for all tested genes and
  returns normalized enrichment scores and leading-edge genes.
- `interpret_enrichment()` uses columns `pathway` and `p_adjust`. Optional
  `genes` contains comma-separated IDs. Map columns before passing results
  from another enrichment package.
- `interpret_ensemble()` uses a strict majority of distinct models by default.
  Shared mechanisms require the same normalized mechanism text and shared
  prompt-grounded evidence. The per-model outputs retain alternative hypotheses.
- The audit reports unsupported pathway names and driver genes against what
  the model actually saw. Evidence grades are model assessments. Enrichment
  and model agreement alone do not establish causal mechanisms.

## Benchmark

A clean run downloads checksum-pinned GEO input files for six perturbations,
fits the sample designs, and generates mixed, up and down ORA baselines.
The optional live benchmark uses a fixed scoring protocol and requires model
credentials. See the [benchmark instructions](benchmark/README.md) and
[report source](benchmark/05_report.Rmd).

## Shiny interface

Install the optional Shiny dependency with `install.packages("shiny")`, then
open the local app with `DEGinsight::run_shiny()`.

Choose the bundled example or upload a CSV/TSV enrichment, differential
expression or signed-ranking table. ORA and GSEA uploads also take a GMT
file. Preranked GSEA requires the optional Bioconductor package `fgsea`.
Use **Prepare evidence** to compute or load enrichment and inspect the prompt.
This step runs locally. **Run LLM** sends the prepared prompt to the
OpenAI-compatible endpoint configured by `LLM_BASE_URL` and `LLM_API_KEY`.
The app presents structured hypotheses and citation checks, with downloads
for enrichment CSV and interpretation JSON.

## License

MIT
