#!/usr/bin/env bash
# Stage report: MultiQC (QC/alignment/RNA) + per-patient self-contained HTML report.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "13_report"

ensure_image "$IMG_FASTQC"; ensure_image "$IMG_PY"
MQ="$RESULTS_DIR/report/multiqc"; mkdir -p "$MQ"
info "MultiQC"
run_logged "$LOG_DIR/multiqc_all.log" \
  dk "$IMG_FASTQC" -m 4g -- multiqc -f -q -o "$MQ" --title "Melanoma pipeline QC" \
    "$RESULTS_DIR/qc" "$RESULTS_DIR/alignment" "$RESULTS_DIR/rna"

for P in $(patients); do
  ctx "$P"; mkdir -p "$D_REP"
  info "$P: building HTML report"
  run_logged "$LOG_DIR/report_${P}.log" \
    dk "$IMG_PY" -m "$PY_DOCKER_MEM" -- python3 "$PROJECT_ROOT/scripts/py/build_report.py" \
      --patient "$P" --results-dir "$RESULTS_DIR" --out "$D_REP/${P}.report.html" \
      --multiqc "$MQ/multiqc_report.html" --aligner "$ALIGNER" --gencode "$GENCODE_RELEASE" --ensembl "$ENSEMBL_RELEASE" \
      --algorithms "$PVAC_ALGORITHMS"
  mark_done "report_${P}"
  info "$P: $D_REP/${P}.report.html"
done
info "done"
