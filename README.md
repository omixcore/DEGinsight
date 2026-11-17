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

## Benchmark

Six public perturbation RNA-seq datasets in which the mechanism is known
independently of the transcriptome (mTORC1, glucocorticoid receptor, HDAC,
BCR-ABL, TNFR1, TLR4), three pre-registered mechanism components per dataset,
and three measured quantities: mechanism recovery, hallucination rate under a
mechanical term audit, and token cost.

The benchmark's load-bearing control is **ORA-only**: the enrichment table the
model is handed is itself scored against the ground truth. Plain ORA places a
canonical pathway in the top five for 4/6 datasets and in the top ten for 5/6
— so a claim that "the model found the pathway" mostly measures the enrichment
step. The interesting quantity is component-level mechanism recovery
(0.00–0.67 across completed runs, $0.005 per dataset median) and the prompt
design that maximises it.

Two results contradict standard practice: few-shot exemplars *reduce* recovery
(0.80 → 0.64) while raising tokens by 31%, and widening the pathway table from
10 to 30 terms is worth more than any instruction we tried.

| document | contents |
|---|---|
| [`benchmark/05_report.Rmd`](benchmark/05_report.Rmd) | full benchmark report; renders to HTML from the CSVs |
| [`paper/preprint.Rmd`](paper/preprint.Rmd) | bioRxiv-style manuscript draft ([docx](paper/preprint.docx) alongside) |
| [`benchmark/README.md`](benchmark/README.md) | how to reproduce every number |

The ORA control (`benchmark/04b_baseline_ora.R`) and both documents render
without an API key. The LLM grid is partially complete; the coverage table at
the top of each document states exactly how many cells are filled, and no
number is imputed.

## Citation

DEGinsight is under active development; a preprint is in preparation.
Until then, please cite the GitHub repository.

## License

MIT
