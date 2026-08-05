#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python3 "$ROOT/run_all.py" "$ROOT/data/otu_table2_with_taxonomy.xlsx" --output-dir "$ROOT/outputs"
