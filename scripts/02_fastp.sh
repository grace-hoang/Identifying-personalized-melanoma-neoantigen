#!/usr/bin/env bash
# Stage QC (2/2): adapter/quality trimming with fastp for WES (tumor, normal) and RNA.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "02_fastp"

ensure_image "$IMG_FASTP"
QC_FASTP="$RESULTS_DIR/qc/fastp"; mkdir -p "$QC_FASTP"

for i in "${!SS_PATIENT[@]}"; do
  P="${SS_PATIENT[i]}"; S="${SS_SAMPLE[i]}"; lane="${SS_LANE[i]}"; fq1="${SS_FQ1[i]}"; fq2="${SS_FQ2[i]}"
  ctx "$P"; mkdir -p "$D_TRIM"
  name="${P}_${S}_L${lane}"
  o1="$D_TRIM/${name}_R1.fastq.gz"; o2="$D_TRIM/${name}_R2.fastq.gz"; js="$QC_FASTP/${name}.fastp.json"
  unit="fastp_${name}"
  if downstream_complete "$P" "$S" && [[ -s "$js" ]]; then info "skip (downstream outputs exist): $name"; continue; fi
  if is_done "$unit" "$o1" "$o2" "$js"; then info "skip (done): $name"; continue; fi
  require_file "$fq1"; require_file "$fq2"
  require_free_gb 25 "$RESULTS_DIR"
  info "fastp $name"
  t1="$D_TRIM/${name}_R1.partial.fastq.gz"; t2="$D_TRIM/${name}_R2.partial.fastq.gz"
  run_logged "$LOG_DIR/fastp_${name}.log" \
    dk "$IMG_FASTP" -m "$FASTP_MEM" -c "$FASTP_THREADS" -- fastp \
      -i "$fq1" -I "$fq2" -o "$t1" -O "$t2" \
      --thread "$FASTP_THREADS" -q "$FASTP_QUALITY" --length_required "$FASTP_LENGTH_REQUIRED" \
      "${FASTP_EXTRA[@]}" --json "$js.partial" --html "$QC_FASTP/${name}.fastp.html" --report_title "$name"
  if (( ! DRY_RUN )); then
    mv -f "$t1" "$o1"; mv -f "$t2" "$o2"; mv -f "$js.partial" "$js"
  fi
  mark_done "$unit"
done
info "done"
