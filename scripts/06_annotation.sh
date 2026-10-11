#!/usr/bin/env bash
# Stage variants (2/2): split multiallelics, left-align, annotate with Ensembl VEP
# (+ pVACtools Frameshift/Wildtype plugins so peptides can be built downstream).
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "06_annotation"

ensure_image "$IMG_BCFTOOLS"; ensure_image "$IMG_VEP"
require_dir "$VEP_CACHE_DIR/homo_sapiens/${ENSEMBL_RELEASE}_GRCh38" "run ./run.sh setup --components vep"
require_file "$VEP_PLUGIN_DIR/Wildtype.pm"; require_file "$VEP_PLUGIN_DIR/Frameshift.pm"

for P in $(patients); do
  ctx "$P"
  if is_done "annotation_${P}" "$VCF_VEP" "$VCF_VEP.tbi"; then info "skip (done): $P"; continue; fi
  require_file "$VCF_PASS" "run stage variants step 05"
  info "normalising $P"
  run_logged "$LOG_DIR/bcftools_norm_${P}.log" \
    dk "$IMG_BCFTOOLS" -m 4g -e "IN=$VCF_PASS" -e "OUT=$VCF_NORM" -e "REF=$REF_FASTA" -- bash -c '
      set -euo pipefail
      bcftools norm -m -any -c w -f "$REF" -Oz -o "$OUT" "$IN"
      bcftools index -t -f "$OUT"
      bcftools stats "$OUT" > "${OUT%.vcf.gz}.stats.txt"'

  info "VEP annotation $P"
  pick=(); (( VEP_PICK )) && pick=(--pick)
  vout="${VCF_VEP%.gz}"
  run_logged "$LOG_DIR/vep_${P}.log" \
    dk "$IMG_VEP" -m "$VEP_DOCKER_MEM" -c "$VEP_FORK" -- vep \
      --input_file "$VCF_NORM" --output_file "$vout" --format vcf --vcf --force_overwrite --no_stats \
      --symbol --terms SO --tsl --mane_select --canonical --biotype --hgvs --fasta "$REF_FASTA" \
      --offline --cache --dir_cache "$VEP_CACHE_DIR" --cache_version "$ENSEMBL_RELEASE" \
      --assembly GRCh38 --species homo_sapiens \
      --plugin Frameshift --plugin Wildtype --dir_plugins "$VEP_PLUGIN_DIR" \
      --fork "$VEP_FORK" --buffer_size 2000 "${pick[@]}"
  run_logged "$LOG_DIR/bcftools_bgzip_${P}.log" \
    dk "$IMG_BCFTOOLS" -m 2g -e "IN=$vout" -e "OUT=$VCF_VEP" -- bash -c '
      set -euo pipefail
      bcftools view -Oz -o "$OUT" "$IN"; bcftools index -t -f "$OUT"'
  discard "$vout"
  mark_done "annotation_${P}"
  info "$P: annotated VCF -> $VCF_VEP"
done
info "done"
