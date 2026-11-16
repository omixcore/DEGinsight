# DEGinsight

**LLM-assisted pathway and mechanism interpretation of bulk RNA-seq differential expression results.**

DEGinsight bridges classical enrichment analysis and large language models:
it grounds LLM reasoning in your actual enrichment statistics, validates
every cited pathway against controlled vocabularies to quantify
hallucinations, and merges multi-model outputs into evidence-graded
mechanistic hypotheses — all reproducible by design (temperature 0 +
response caching).

## Why not just paste results into ChatGPT?

| | Manual copy-paste | DEGinsight |
|---|---|---|
| Prompt construction | ad hoc, lossy | ontology-constrained, includes stats + driver genes |
| Hallucination control | none | every cited term validated against GO/KEGG/input vocabulary |
| Reproducibility | none | temperature 0, on-disk response cache, versioned prompts |
| Multi-model evidence | single opinion | elite-ensemble majority vote with agreement metrics |
| Cost control | unknown | `estimate_cost()` before you spend |

## Installation

```r
# install.packages("remotes")
remotes::install_github("omixcore/DEGinsight")
```

## 5-minute example

```r
library(DEGinsight)
Sys.setenv(LLM_API_KEY = "sk-...")                       # any OpenAI-compatible key
Sys.setenv(LLM_BASE_URL = "https://api.deepseek.com/v1") # or OpenAI / Moonshot / local Ollama

data(example_data)
res <- interpret_enrichment(
  example_ora,
  contrast = "TNF-stimulated vs control HeLa cells",
  species  = "Homo sapiens",
  model    = "deepseek-chat"
)
res                                   # hypotheses + hallucination audit

# ensemble mode: majority vote across frontier models
out <- interpret_deg(example_deg, gene_col = "symbol", padj_col = "padj",
                     logfc_col = "log2FoldChange",
                     gene_sets = example_gene_sets,
                     contrast = "TNF vs control", species = "Homo sapiens",
                     models = c("deepseek-chat", "kimi-k2", "gpt-4o-mini"))

render_report(out$ora, out$interpretation, "report.html")
```

## Features

- **Unified backend** — one function for OpenAI / DeepSeek / Moonshot /
  Azure / local Ollama; swap models by changing a string.
- **Ontology-grounded prompting** — the model is shown the exact enriched
  terms and forbidden from inventing new ones.
- **Hallucination audit** — `validate_terms()` returns per-hypothesis and
  global hallucination rates (a benchmark headline metric).
- **Ensemble interpretation** — pathway-level majority voting with
  Jaccard agreement matrices across models.
- **Cost estimation** — know the price of a benchmark run beforehand.
- **Lightweight** — core ORA has zero Bioconductor dependencies; works with
  `clusterProfiler` output when you need full databases.
- **Privacy option** — point `LLM_BASE_URL` at a local Ollama instance for
  controlled data.

## Benchmark (in progress)

We are benchmarking 6+ LLMs on mechanistic recovery using perturbation
datasets with known ground-truth mechanisms (LINCS, GEO knockdown series).
Metrics: mechanism recall, hallucination rate, confidence calibration and
cost-performance trade-off. Results will be posted on bioRxiv.

## Citation

DEGinsight is under active development; a preprint is in preparation.
Until then, please cite the GitHub repository.

## License

MIT
