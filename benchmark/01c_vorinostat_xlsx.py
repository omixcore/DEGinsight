"""01c_vorinostat_xlsx.py — convert GSE319257 per-sample xlsx to one counts CSV.

GEO deposited this series as one xlsx per sample (no combined matrix), with a
single sheet holding 45k genes. Column A = gene symbol, column D =
"Expression value" (read count).

openpyxl needs minutes per file on these workbooks, so this script streams
the sheet XML directly with ElementTree, which takes a few seconds.

Usage:  python benchmark/01c_vorinostat_xlsx.py
Input:  benchmark/geo_cache/vorinostat_xlsx/*.xlsx
Output: benchmark/geo_cache/GSE319257_counts.csv
"""
import csv
import glob
import os
import re
import sys
import zipfile
from xml.etree.ElementTree import iterparse

SRC = os.path.join("benchmark", "geo_cache", "vorinostat_xlsx")
OUT = os.path.join("benchmark", "geo_cache", "GSE319257_counts.csv")
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def shared_strings(zf):
    out = []
    with zf.open("xl/sharedStrings.xml") as fh:
        for _, el in iterparse(fh, events=("end",)):
            if el.tag == NS + "si":
                out.append("".join(t.text or "" for t in el.iter(NS + "t")))
                el.clear()
    return out


def sheet_path(zf):
    names = [n for n in zf.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", n)]
    if not names:
        sys.exit("no worksheet in archive")
    return sorted(names)[0]


def read_sample(path):
    """Return {gene: count} using column A (symbol) and column D (count)."""
    with zipfile.ZipFile(path) as zf:
        try:
            sst = shared_strings(zf)
        except KeyError:
            sst = []
        genes, row = {}, 0
        with zf.open(sheet_path(zf)) as fh:
            for _, el in iterparse(fh, events=("end",)):
                if el.tag == NS + "row":
                    row += 1
                    a = d = None
                    for c in el.iter(NS + "c"):
                        ref = c.get("r", "")
                        col = re.sub(r"\d", "", ref)
                        if col not in ("A", "D"):
                            continue
                        kind = c.get("t")
                        if kind == "inlineStr":
                            node = c.find(NS + "is")
                            val = "".join(t.text or "" for t in node.iter(NS + "t")) \
                                if node is not None else None
                        else:
                            v = c.find(NS + "v")
                            if v is None or v.text is None:
                                continue
                            val = sst[int(v.text)] if kind == "s" else v.text
                        if val is None:
                            continue
                        if col == "A":
                            a = val
                        else:
                            d = val
                    if a is not None and d is not None and row > 1:
                        try:
                            genes[a] = float(d)
                        except ValueError:
                            pass
                    el.clear()
    return genes


files = sorted(glob.glob(os.path.join(SRC, "*.xlsx")))
if not files:
    sys.exit(f"no xlsx found in {SRC}")

samples, series = [], {}
for path in files:
    name = re.sub(r"^GSM\d+_", "", os.path.basename(path)).replace(".xlsx", "")
    series[name] = read_sample(path)
    samples.append(name)
    print(f"{name}: {len(series[name])} genes", flush=True)

genes = sorted(set().union(*[set(s) for s in series.values()]))
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene"] + samples)
    for g in genes:
        w.writerow([g] + [series[s].get(g, 0.0) for s in samples])

print(f"wrote {OUT}: {len(genes)} genes x {len(samples)} samples")
