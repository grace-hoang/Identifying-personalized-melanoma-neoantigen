#!/usr/bin/env bash
# Stage neoantigen (2/3): composite ranking of pVACseq candidates (binding, agretopicity,
# clonality, expression, RNA support, HLA promiscuity) and selection for construct design.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "11_ranking"

ensure_image "$IMG_PY"
for P in $(patients); do
  ctx "$P"
  if is_done "ranking_${P}" "$RANKED" "$DESIGN_INPUT"; then info "skip (done): $P"; continue; fi
  allep="$(ls -1 "$PVAC_OUT"/MHC_Class_I/*.all_epitopes.tsv 2>/dev/null | head -n 1 || true)"
  if [[ -z "$allep" ]]; then
    if (( DRY_RUN )); then allep="$PVAC_OUT/MHC_Class_I/${SM_TUMOR}.all_epitopes.tsv"; else die "$P: no pVACseq all_epitopes.tsv — run step 10"; fi
  fi
  info "$P: ranking $(basename "$allep")"
  run_logged "$LOG_DIR/ranking_${P}.log" \
    dk "$IMG_PY" -m "$PY_DOCKER_MEM" -- python3 "$PROJECT_ROOT/scripts/py/rank_neoantigens.py" \
      --all-epitopes "$allep" --patient "$P" --out-dir "$D_NEO" \
      --ic50-max "$RANK_IC50_MAX" --percentile-max "$RANK_PERCENTILE_MAX" \
      --min-fold-change "$RANK_MIN_FOLD_CHANGE" --min-tumor-vaf "$RANK_MIN_TUMOR_VAF" \
      --min-tumor-dp "$RANK_MIN_TUMOR_DP" --max-normal-vaf "$RANK_MAX_NORMAL_VAF" \
      --min-expr-tpm "$RANK_MIN_EXPR_TPM" --top-n "$RANK_TOP_N" --select-n "$DESIGN_N_EPITOPES" \
      ${TUMOR_PURITY:+--purity "$TUMOR_PURITY"}
  mark_done "ranking_${P}"
done
info "done"
