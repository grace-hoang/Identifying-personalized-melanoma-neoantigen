#!/usr/bin/env bash
# Stage neoantigen (3/3): in-silico design of a candidate multi-epitope mRNA construct
# (signal peptide - epitopes/linkers - MITD, codon optimisation, UTRs, poly(A), QC, RNAfold).
# RESEARCH USE ONLY — computational candidate, not a validated or clinical product.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "12_mrna_design"

ensure_image "$IMG_PY"
for P in $(patients); do
  ctx "$P"; mkdir -p "$D_MRNA"
  if is_done "mrna_${P}" "$D_MRNA/${P}.construct.mrna.fa" "$D_MRNA/${P}.design_report.json"; then info "skip (done): $P"; continue; fi
  require_file "$DESIGN_INPUT" "run step 11_ranking"
  info "$P: designing construct"
  run_logged "$LOG_DIR/mrna_design_${P}.log" \
    dk "$IMG_PY" -m "$PY_DOCKER_MEM" -- python3 "$PROJECT_ROOT/scripts/py/design_mrna.py" \
      --patient "$P" --design-input "$DESIGN_INPUT" --vep-vcf "$VCF_VEP" --design-dir "$DESIGN_DIR" \
      --out-dir "$D_MRNA" --max-epitopes "$DESIGN_N_EPITOPES" --flank "$DESIGN_FLANK" \
      --linker-short "$DESIGN_LINKER_SHORT" --linker-long "$DESIGN_LINKER_LONG" \
      --signal "$DESIGN_INCLUDE_SIGNAL" --mitd "$DESIGN_INCLUDE_MITD" --padre "$DESIGN_INCLUDE_PADRE" \
      --polya-len "$DESIGN_POLYA_LEN" --seed "$DESIGN_SEED" --iters "$DESIGN_ANNEAL_ITERS" \
      --gc-min "$DESIGN_GC_MIN" --gc-max "$DESIGN_GC_MAX" --avoid-sites "$DESIGN_AVOID_SITES" \
      --t7-promoter "$DESIGN_T7_PROMOTER"
  mark_done "mrna_${P}"
  info "$P: construct written to $D_MRNA"
done
info "done"
