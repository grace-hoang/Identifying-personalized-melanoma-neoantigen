#!/usr/bin/env bash
# Stage QC (1/2): FastQC on every raw FASTQ + MultiQC summary of raw reads.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "01_fastqc"

QC_RAW="$RESULTS_DIR/qc/fastqc_raw"
mkdir -p "$QC_RAW"
ensure_image "$IMG_FASTQC"

todo=(); todo_units=()
for i in "${!SS_PATIENT[@]}"; do
  for f in "${SS_FQ1[i]}" "${SS_FQ2[i]}"; do
    require_file "$f"
    base="$(fq_base "$f")"
    if is_done "fastqc_raw_${base}" "$QC_RAW/${base}_fastqc.zip"; then
      info "skip (done): FastQC $base"
    else
      todo+=("$f"); todo_units+=("fastqc_raw_${base}")
    fi
  done
done

if (( ${#todo[@]} > 0 )); then
  info "running FastQC on ${#todo[@]} file(s)"
  run_logged "$LOG_DIR/fastqc_raw.log" \
    dk "$IMG_FASTQC" -m "$FASTQC_MEM" -c "$FASTQC_THREADS" -- \
       fastqc -t "$FASTQC_THREADS" --memory 1024 -o "$QC_RAW" -d "$TMP_DIR" "${todo[@]}"
  for u in "${todo_units[@]}"; do mark_done "$u"; done
fi

if is_done "multiqc_raw" "$RESULTS_DIR/qc/multiqc_raw/multiqc_report.html" && (( ${#todo[@]} == 0 )); then
  info "skip (done): MultiQC raw"
else
  run_logged "$LOG_DIR/multiqc_raw.log" \
    dk "$IMG_FASTQC" -m 4g -- multiqc -f -q -o "$RESULTS_DIR/qc/multiqc_raw" "$QC_RAW"
  mark_done "multiqc_raw"
fi
info "done"
