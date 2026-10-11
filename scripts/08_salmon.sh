#!/usr/bin/env bash
# Stage RNA (2/2): Salmon transcript quantification -> transcript- and gene-level TPM tables.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "08_salmon"

ensure_image "$IMG_SALMON"; ensure_image "$IMG_PY"
require_file "$SALMON_INDEX_DIR/versionInfo.json" "run ./run.sh setup --components salmon"
require_file "$TX2GENE"

for P in $(patients); do
  ctx "$P"
  has_sample "$P" RNA || { warn "$P: no RNA row, skipping"; continue; }
  mkdir -p "$D_RNA/salmon" "$D_RNA/expression"
  TXT="$D_RNA/expression/${P}.transcript_expression.tsv"; GXT="$D_RNA/expression/${P}.gene_expression.tsv"
  if is_done "salmon_${P}" "$QUANT_SF" "$TXT" "$GXT"; then info "skip (done): $P"; else
    r1=(); r2=()
    while IFS=$'\t' read -r lane _ _; do
      a="$D_TRIM/${P}_RNA_L${lane}_R1.fastq.gz"; b="$D_TRIM/${P}_RNA_L${lane}_R2.fastq.gz"
      require_file "$a" "run stage qc first"; require_file "$b"
      r1+=("$a"); r2+=("$b")
    done < <(rows "$P" RNA)
    info "salmon quant $P"
    run_logged "$LOG_DIR/salmon_${P}.log" \
      dk "$IMG_SALMON" -m "$SALMON_DOCKER_MEM" -c "$SALMON_THREADS" -- salmon quant \
        -i "$SALMON_INDEX_DIR" -l A -1 "${r1[@]}" -2 "${r2[@]}" -p "$SALMON_THREADS" \
        --validateMappings --seqBias --gcBias -o "$D_RNA/salmon"
    run_logged "$LOG_DIR/expr_tables_${P}.log" \
      dk "$IMG_PY" -m "$PY_DOCKER_MEM" -- python3 "$PROJECT_ROOT/scripts/py/expr_tables.py" \
        --quant "$QUANT_SF" --tx2gene "$TX2GENE" --out-tx "$TXT" --out-gene "$GXT"
    mark_done "salmon_${P}"
  fi
  if (( CLEAN_TRIMMED_FASTQ )) && is_done "star_${P}" "$BAM_RNA" && is_done "salmon_${P}" "$QUANT_SF"; then
    discard "$D_TRIM"/${P}_RNA_L*_R[12].fastq.gz
  fi
done
info "done"
