#!/usr/bin/env bash
# 99_finish.sh — resume the benchmark after topping up API credit.
#
# Fills the remaining dataset x model cells, adds the cross-model judge score,
# recomputes the ORA control, and re-renders both documents. Everything is
# idempotent: completed cells are skipped, cached responses are reused, so
# re-running this script costs nothing.
#
#   export LLM_API_KEY=sk-...
#   bash benchmark/99_finish.sh
#
# Optional:
#   BENCH_MODELS=deepseek-flash,deepseek-v4-pro   # restrict the grid
#   ABLATION_MODEL=deepseek-v4-pro                # redo the ablation on another model

set -euo pipefail
cd "$(dirname "$0")/.."

# --- locate R -------------------------------------------------------------
if ! command -v Rscript >/dev/null 2>&1; then
  for d in "/c/Program Files/R"/R-*/bin/x64; do
    [ -x "$d/Rscript.exe" ] && export PATH="$d:$PATH" && break
  done
fi
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found on PATH"; exit 1; }

# pandoc ships with Positron/RStudio on this machine; rmarkdown needs to know
if [ -z "${RSTUDIO_PANDOC:-}" ]; then
  for p in "/c/Program Files/Positron/bin/pandoc" \
           "/c/Program Files/RStudio/bin/pandoc"; do
    [ -d "$p" ] && export RSTUDIO_PANDOC="$p" && break
  done
fi
export R_LIBS_USER="${R_LIBS_USER:-G:/WorkBuddy/2026-11-16-20-16-33/Rlib}"

# --- precondition ---------------------------------------------------------
if [ -z "${LLM_API_KEY:-}" ]; then
  echo "LLM_API_KEY is not set — export it first (see benchmark/models.csv)."
  exit 1
fi

probe=$(curl -s -m 20 https://api.deepseek.com/user/balance \
          -H "Authorization: Bearer $LLM_API_KEY" || true)
case "$probe" in
  *'"is_available":true'*) echo "balance OK" ;;
  *'"is_available":false'*) echo "API reports no available balance — top up and retry:"; echo "$probe"; exit 1 ;;
  *) echo "could not read balance (endpoint may not be DeepSeek); continuing anyway" ;;
esac

# --- pipeline -------------------------------------------------------------
echo "== 1/5 remaining dataset x model cells =="
Rscript benchmark/02_run_benchmark.R

echo "== 2/5 scoring + cross-model judge =="
RUN_JUDGE=1 JUDGE_MODE=cross Rscript benchmark/04_score.R

echo "== 3/5 ORA-only control (free) =="
Rscript benchmark/04b_baseline_ora.R

echo "== 4/5 documents =="
( cd benchmark && Rscript -e 'rmarkdown::render("05_report.Rmd", output_file="05_report.html", quiet=TRUE)' )
( cd paper     && Rscript -e 'rmarkdown::render("preprint.Rmd",  output_file="preprint.html",  quiet=TRUE)' \
              && Rscript -e 'rmarkdown::render("preprint.Rmd",  output_file="preprint.docx",  quiet=TRUE)' )

echo "== 5/5 done =="
echo "Reports:  benchmark/05_report.html"
echo "Manuscript: paper/preprint.html  paper/preprint.docx"
echo "Tables:   benchmark/scores/{benchmark_summary,model_ranking,pathway_recall}.csv"
