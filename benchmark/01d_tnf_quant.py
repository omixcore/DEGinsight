"""01d_tnf_quant.py — Salmon quant.sf (GSE305796) -> gene-level count matrix.

The series deposits one Salmon `quant.sf.gz` per sample (transcript-level
NumReads, GENCODE transcript IDs), so gene symbols must be recovered before
ORA can be run. Transcript -> symbol mapping is fetched from mygene.info
(no API key, batch queries), and counts are summed per symbol.

Samples used (Untreated vs 30 ng/mL TNF-alpha, three donor lines):
  GSM9185037 C156 / GSM9185038 T156
  GSM9185044 C67  / GSM9185045 T67
  GSM9185051 C840 / GSM9185052 T840

Usage:  python benchmark/01d_tnf_quant.py
Input:  benchmark/geo_cache/tnf/*_quant.sf.gz
Output: benchmark/geo_cache/GSE305796_counts.csv
"""
import glob
import gzip
import json
import os
import sys
import urllib.request

SRC = os.path.join("benchmark", "geo_cache", "tnf")
OUT = os.path.join("benchmark", "geo_cache", "GSE305796_counts.csv")

WANTED = [
    ("C156", "GSM9185037_C156_quant.sf.gz"),
    ("T156", "GSM9185038_T156_quant.sf.gz"),
    ("C67",  "GSM9185044_C67_quant.sf.gz"),
    ("T67",  "GSM9185045_T67_quant.sf.gz"),
    ("C840", "GSM9185051_C840_quant.sf.gz"),
    ("T840", "GSM9185052_T840_quant.sf.gz"),
]

if not all(os.path.exists(os.path.join(SRC, f)) for _, f in WANTED):
    have = sorted(glob.glob(os.path.join(SRC, "*.gz")))
    sys.exit(f"missing quant files in {SRC}; found {len(have)}")

# ---- read NumReads per transcript -----------------------------------------
samples, series = [], {}
for label, fname in WANTED:
    counts = {}
    with gzip.open(os.path.join(SRC, fname), "rt") as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            counts[parts[0].split(".")[0]] = float(parts[4])
    series[label] = counts
    samples.append(label)
    print(f"{label}: {len(counts)} transcripts", flush=True)

transcripts = sorted(set().union(*[set(s) for s in series.values()]))
# drop transcripts that are essentially unexpressed in every sample
transcripts = [t for t in transcripts if sum(series[s].get(t, 0.0) for s in samples) >= 10]
print(f"{len(transcripts)} transcripts after expression filter", flush=True)

# ---- transcript -> symbol via mygene.info --------------------------------
def query_batch(batch):
    # mygene's batch POST takes the terms in `q` (not `qterms`)
    payload = json.dumps({
        "q": batch,
        "scopes": "ensembl.transcript",
        "fields": "symbol",
        "species": "human",
    }).encode()
    req = urllib.request.Request(
        "https://mygene.info/v3/query",
        data=payload,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.load(fh)


mapping = {}
B = 1000
for i in range(0, len(transcripts), B):
    batch = transcripts[i:i + B]
    try:
        res = query_batch(batch)
    except Exception as exc:  # noqa: BLE001
        print(f"  batch {i} failed: {exc}", flush=True)
        continue
    for t, hit in zip(batch, res):
        if isinstance(hit, dict) and hit.get("symbol") and not hit.get("notfound"):
            mapping[t] = hit["symbol"]
    print(f"  mapped {len(mapping)}/{min(i + B, len(transcripts))}", flush=True)

print(f"mapped {len(mapping)} / {len(transcripts)} transcripts to symbols")

# ---- aggregate counts per symbol -----------------------------------------
genes = {}
for t in transcripts:
    sym = mapping.get(t)
    if not sym:
        continue
    row = genes.setdefault(sym, [0.0] * len(samples))
    for j, s in enumerate(samples):
        row[j] += series[s].get(t, 0.0)

with open(OUT, "w", newline="") as fh:
    fh.write("gene," + ",".join(samples) + "\n")
    for sym in sorted(genes):
        fh.write(sym + "," + ",".join(f"{v:g}" for v in genes[sym]) + "\n")

print(f"wrote {OUT}: {len(genes)} genes x {len(samples)} samples")
