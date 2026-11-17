#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p geo_cache
while read -r sha url name; do
  file="geo_cache/$name"
  if printf '%s  %s\n' "$sha" "$file" | sha256sum --check --status; then
    continue
  fi
  curl --fail --silent --show-error --location --retry 3 --max-time 600 "$url" -o "$file.part"
  printf '%s  %s\n' "$sha" "$file.part" | sha256sum --check
  mv "$file.part" "$file"
done < sources.tsv
