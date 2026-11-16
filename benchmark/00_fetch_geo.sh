#!/usr/bin/env bash
# 00_fetch_geo.sh — fetch raw supplementary files for the benchmark datasets.
#
# NCBI throttles each TCP connection to ~15 KB/s, so downloading a 139 MB
# archive over a single connection takes hours. This script splits every
# archive into PARTS byte ranges and fetches them concurrently, then
# concatenates. Each part resumes from where it stopped, so the script is
# safe to interrupt and re-run.
#
# Usage:  bash benchmark/00_fetch_geo.sh
# Output: benchmark/geo_cache/

set -u
cd "$(dirname "$0")/geo_cache" || exit 1
PARTS=${PARTS:-10}

dl_part () {  # dl_part <url> <start> <end> <target>
  local url="$1" s="$2" e="$3" tgt="$4" sz need try
  need=$(( e - s + 1 ))
  for try in $(seq 1 40); do
    sz=$(stat -c%s "$tgt" 2>/dev/null || echo 0)
    [ "$sz" -ge "$need" ] && return 0
    curl -sfk --max-time 300 -r $(( s + sz ))-"$e" "$url" >> "$tgt" || true
  done
  sz=$(stat -c%s "$tgt" 2>/dev/null || echo 0)
  [ "$sz" -ge "$need" ]
}

fetch_par () {  # fetch_par <url> <expected_bytes> <outfile>
  local url="$1" want="$2" out="$3" chunk i s e sz
  sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [ "$sz" -eq "$want" ]; then echo "OK       $out ($sz bytes)"; return 0; fi

  rm -f "$out"
  chunk=$(( (want + PARTS - 1) / PARTS ))
  for i in $(seq 0 $(( PARTS - 1 ))); do
    s=$(( i * chunk )); e=$(( s + chunk - 1 ))
    [ "$e" -ge "$want" ] && e=$(( want - 1 ))
    dl_part "$url" "$s" "$e" "$out.part$i" &
  done
  wait

  : > "$out"
  for i in $(seq 0 $(( PARTS - 1 ))); do cat "$out.part$i" >> "$out" 2>/dev/null; done
  rm -f "$out".part*
  sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [ "$sz" -eq "$want" ]; then echo "OK       $out ($sz bytes)"; else
    echo "FAILED   $out got $sz / want $want"; return 1; fi
}

base="https://ftp.ncbi.nlm.nih.gov/geo/series"

# sizes verified against Content-Length on 2026-11-16
fetch_par "$base/GSE318nnn/GSE318529/suppl/GSE318529_Raw_Counts.csv.gz"       1321965 GSE318529_Raw_Counts.csv.gz
fetch_par "$base/GSE304nnn/GSE304202/suppl/GSE304202_normalized_counts.csv.gz" 1665726 GSE304202_normalized_counts.csv.gz
fetch_par "$base/GSE297nnn/GSE297509/suppl/GSE297509_RAW.tar"                 28958720 GSE297509_RAW.tar
fetch_par "$base/GSE305nnn/GSE305796/suppl/GSE305796_RAW.tar"                 49848320 GSE305796_RAW.tar
fetch_par "$base/GSE319nnn/GSE319257/suppl/GSE319257_RAW.tar"                139038720 GSE319257_RAW.tar

echo "--- verifying archives ---"
for f in *.tar; do
  [ -e "$f" ] || continue
  if tar -tf "$f" >/dev/null 2>&1; then echo "OK   $f"; else echo "BAD  $f"; fi
done
