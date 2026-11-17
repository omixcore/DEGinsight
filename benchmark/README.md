# DEGinsight benchmark

The benchmark measures recovery of mechanisms associated with six known
perturbations. Models see the enrichment table, species, a generic
`treated vs control` contrast and, where enabled, tissue or cell line. The
known perturbagen is omitted from the prompt. The few-shot example uses
fictional terms so it does not supply an expected benchmark mechanism.

## Reproduce from the repository root

Requirements include R with limma, edgeR and digest, Python 3.9+ with
openpyxl, Bash, curl and sha256sum. Pandoc and rmarkdown are used for the report.

```bash
Rscript benchmark/run_all.R
Rscript -e 'rmarkdown::render("benchmark/05_report.Rmd")'
```

Inputs and SHA-256 hashes are in [sources.tsv](sources.tsv). The pipeline
writes DEG tables with statistical results, sample manifests, designs and R
session information to `benchmark/deg_tables/`. The deterministic ORA baseline
is written to `benchmark/revised/scores/baseline_ora.csv`.

| Dataset | Expression input | Model |
| --- | --- | --- |
| [GSE318529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318529) Rapamycin | Raw counts, Stage 4 15-C | TMM + voom, paired by donor (3 pairs) |
| [GSE304202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304202) Dexamethasone | Submitted normalized counts, log2(x + 1) | limma-trend, paired by published matched design and submitted IDs (6 pairs) |
| [GSE321679](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE321679) Imatinib | Submitted voom-scale logCPM, K562 WT | limma-trend, GEO batch blocks (3 pairs). Original voom weights unavailable |
| [GSE319257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319257) Vorinostat | Unique reads from the mock-transfected A498 arm | TMM + voom (6 per group). Supplement gene IDs are human ENSG IDs |
| [GSE305796](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305796) TNF | Salmon transcript TPM summarized using Ensembl 97 and scaled to library size | TMM + voom, paired by donor (3 pairs) |
| [GSE297509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297509) LPS | Unique reads from the control-siRNA arm | TMM + voom, paired by donor (2 pairs) |

Count matrices are filtered with edgeR `filterByExpr()` before normalization.
Normalized Dexamethasone counts retain genes with at least one normalized
count in at least six samples. DE is fit by gene ID and adjusted by BH.
ORA maps tested genes to symbols and resolves duplicate symbols by highest
mean expression, with ties resolved by gene ID. It applies adjusted
P < 0.05 and absolute log2 fold change > 1 throughout. Mixed, up and down ORA
each use the same tested-symbol universe and their own BH correction. LPS has
one residual degree of freedom, which limits precision.

## Live benchmark and scoring

Configure the credentials and endpoints in [models.csv](models.csv) through
the listed environment variables. `BENCH_MODELS` defaults to
`deepseek-flash,deepseek-v4-pro`. `ABLATION_MODEL` defaults to `deepseek-flash`.

```bash
Rscript benchmark/run_all.R --live
# Optional fixed judge. Supply credentials for this distinct model.
RUN_JUDGE=1 JUDGE_MODEL=gpt-4o-mini Rscript benchmark/04_score.R
```

- The main grid evaluates six data sets with each model using up to 20 mixed
  ORA pathways with positive overlaps. The ablation evaluates three data sets
  at 5, 10, 20 and 30 pathways, with and without few-shot and tissue context.
- Pathway recovery compares grounded citations to the ORA prefix with the same
  number of terms. Prompt coverage gives the number of expected pathways
  available to the model. Keyword matches are descriptive, with a separate
  optional fixed-model semantic score. The terminology audit covers named
  pathways and driver genes visible in the prompt.
- Results are stored with a fingerprint of code, configuration, DEG input
  and references. They can be resumed without loading legacy responses.
- The revised scoring reference is in [ground_truth.csv](ground_truth.csv).
  Its ancestry includes edits after the first exploratory model calls. The
  reused studies and their reference terms constitute a retrospective
  benchmark. Independent validation needs a separate held-out dataset set.

The original `results/`, `ablation/` and `scores/` directories and the
[preprint draft](../paper/preprint.Rmd) contain exploratory outputs under the
earlier protocol and are archived for provenance.
