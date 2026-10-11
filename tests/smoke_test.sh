#!/usr/bin/env bash
# Offline smoke test (no Docker, no references): synthetic data through ranking, design and report.
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
python3 "$ROOT/tests/make_synthetic.py" "$T/in" >/dev/null
python3 "$ROOT/scripts/py/rank_neoantigens.py" --all-epitopes "$T/in/syn.all_epitopes.tsv" --patient SYN --out-dir "$T/res/neoantigen/SYN" >/dev/null
python3 "$ROOT/scripts/py/design_mrna.py" --patient SYN --design-input "$T/res/neoantigen/SYN/SYN.design_input.tsv" \
  --vep-vcf "$T/in/syn.vep.vcf.gz" --design-dir "$ROOT/config/design" --out-dir "$T/res/mrna/SYN" --max-epitopes 8 --iters 1500 >/dev/null
python3 "$ROOT/scripts/py/build_report.py" --patient SYN --results-dir "$T/res" --out "$T/res/report/SYN.html" >/dev/null
test -s "$T/res/report/SYN.html" && test -s "$T/res/mrna/SYN/SYN.construct.mrna.fa" && echo "smoke test OK"
