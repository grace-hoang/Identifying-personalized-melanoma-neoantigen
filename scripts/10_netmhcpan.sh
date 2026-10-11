#!/usr/bin/env bash
# Stage neoantigen (1/3): build the pVACseq input (VEP VCF + expression + RNA read counts)
# and run pVACseq (NetMHCpan / NetMHCpanEL / MHCflurry ... as set in PVAC_ALGORITHMS).
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "10_netmhcpan"

ensure_image "$IMG_PVAC"; ensure_image "$IMG_BCFTOOLS"; ensure_image "$IMG_BAMRC"

for P in $(patients); do
  ctx "$P"; mkdir -p "$D_NEO"
  if is_done "pvacseq_${P}" "$PVAC_INPUT" ; then
    if compgen -G "$PVAC_OUT/MHC_Class_I/*all_epitopes.tsv" >/dev/null; then info "skip (done): $P"; continue; fi
  fi
  require_file "$VCF_VEP" "run stage variants (06_annotation)"
  require_file "$HLA_FILE" "run stage hla (09_hla)"
  hla="$(tr -d '[:space:]' < "$HLA_FILE" 2>/dev/null || echo "HLA-A*02:01")"
  info "$P: HLA alleles for prediction: $hla"

  # ---------------- 1. prepare pVACseq input ----------------------------------
  if ! is_done "pvac_input_${P}" "$PVAC_INPUT" "$PVAC_INPUT.tbi"; then
    W="$D_NEO/prep"; mkdir -p "$W"; register_cleanup "$W"
    cur="$W/step0.vcf"
    run_logged "$LOG_DIR/pvac_prep_${P}.log" \
      dk "$IMG_BCFTOOLS" -m 2g -e "IN=$VCF_VEP" -e "OUT=$cur" -- bash -c 'bcftools view -Ov -o "$OUT" "$IN"'

    # expression (transcript + gene TPM from salmon)
    TXT="$D_RNA/expression/${P}.transcript_expression.tsv"; GXT="$D_RNA/expression/${P}.gene_expression.tsv"
    if [[ -s "$TXT" && -s "$GXT" ]]; then
      info "$P: annotating transcript/gene expression"
      run_logged "$LOG_DIR/pvac_prep_${P}.log" \
        dk "$IMG_PVAC" -m 4g -- vcf-expression-annotator "$cur" "$TXT" custom transcript \
          -s "$SM_TUMOR" --id-column transcript_id --expression-column tpm --ignore-ensembl-id-version -o "$W/step1.vcf"
      run_logged "$LOG_DIR/pvac_prep_${P}.log" \
        dk "$IMG_PVAC" -m 4g -- vcf-expression-annotator "$W/step1.vcf" "$GXT" custom gene \
          -s "$SM_TUMOR" --id-column gene_id --expression-column tpm --ignore-ensembl-id-version -o "$W/step2.vcf"
      cur="$W/step2.vcf"
    else
      warn "$P: expression tables missing (run stage rna) — continuing WITHOUT expression annotation"
    fi

    # RNA allele read counts (bam-readcount) -> RAD/RDP/RAF
    if [[ -s "$BAM_RNA" ]]; then
      info "$P: RNA read counts at variant sites"
      rc="$W/rc"; mkdir -p "$rc"
      run_logged "$LOG_DIR/pvac_prep_${P}.log" \
        dk "$IMG_BAMRC" -m 8g -- /usr/bin/python /usr/bin/bam_readcount_helper.py \
          "$cur" "$SM_TUMOR" "$REF_FASTA" "$BAM_RNA" NOPREFIX "$rc"
      snv="$rc/${SM_TUMOR}_bam_readcount_snv.tsv"; ind="$rc/${SM_TUMOR}_bam_readcount_indel.tsv"
      if [[ -s "$snv" ]]; then
        run_logged "$LOG_DIR/pvac_prep_${P}.log" \
          dk "$IMG_PVAC" -m 4g -- vcf-readcount-annotator "$cur" "$snv" RNA -s "$SM_TUMOR" -t snv -o "$W/rc_snv.vcf"
        cur="$W/rc_snv.vcf"
      fi
      if [[ -s "$ind" ]]; then
        run_logged "$LOG_DIR/pvac_prep_${P}.log" \
          dk "$IMG_PVAC" -m 4g -- vcf-readcount-annotator "$cur" "$ind" RNA -s "$SM_TUMOR" -t indel -o "$W/rc_indel.vcf"
        cur="$W/rc_indel.vcf"
      fi
    else
      warn "$P: no RNA BAM — continuing WITHOUT RNA read counts"
    fi

    run_logged "$LOG_DIR/pvac_prep_${P}.log" \
      dk "$IMG_BCFTOOLS" -m 2g -e "IN=$cur" -e "OUT=$PVAC_INPUT" -- bash -c \
        'set -euo pipefail; bcftools view -Oz -o "$OUT" "$IN"; bcftools index -t -f "$OUT"'
    mark_done "pvac_input_${P}"
  fi

  # ---------------- 2. pVACseq ------------------------------------------------
  rm -rf "$PVAC_OUT"; mkdir -p "$PVAC_OUT"
  read -ra algs <<< "$PVAC_ALGORITHMS"
  iedb=()
  mounts=()
  if [[ -n "$PVAC_LOCAL_IEDB_HOST_DIR" ]]; then mounts=(-v "$PVAC_LOCAL_IEDB_HOST_DIR:/iedb:ro"); iedb=(--iedb-install-directory /iedb)
  elif [[ "$PVAC_IEDB_DIR" == "auto" ]]; then
    if (( ! DRY_RUN )) && docker run --rm "$IMG_PVAC" test -d /opt/iedb 2>/dev/null; then iedb=(--iedb-install-directory /opt/iedb)
    else warn "no local IEDB found in the image: pVACseq will call the IEDB REST API (internet needed, slow)"; fi
  elif [[ "$PVAC_IEDB_DIR" != "none" ]]; then iedb=(--iedb-install-directory "$PVAC_IEDB_DIR"); fi
  purity=(); [[ -n "$TUMOR_PURITY" ]] && purity=(--tumor-purity "$TUMOR_PURITY")

  info "$P: pVACseq (${algs[*]}; lengths $PVAC_CLASS_I_LENGTHS)"
  run_logged "$LOG_DIR/pvacseq_${P}.log" \
    dk "$IMG_PVAC" -m "$PVAC_DOCKER_MEM" -c "$PVAC_THREADS" "${mounts[@]}" -- pvacseq run \
      "$PVAC_INPUT" "$SM_TUMOR" "$hla" "${algs[@]}" "$PVAC_OUT" \
      -e1 "$PVAC_CLASS_I_LENGTHS" -t "$PVAC_THREADS" --normal-sample-name "$SM_NORMAL" \
      -b "$PVAC_BINDING_THRESHOLD" --binding-percentile-threshold "$PVAC_PERCENTILE_THRESHOLD" \
      -c "$PVAC_MIN_FOLD_CHANGE" --normal-vaf "$PVAC_NORMAL_VAF_MAX" --tdna-vaf "$PVAC_TUMOR_DNA_VAF_MIN" \
      --expn-val "$PVAC_EXPN_VAL" --pass-only -a sample_name \
      "${iedb[@]}" "${purity[@]}" "${PVAC_EXTRA_ARGS[@]}"
  if (( ! DRY_RUN )) && ! compgen -G "$PVAC_OUT/MHC_Class_I/*all_epitopes.tsv" >/dev/null; then
    die "pVACseq finished but no all_epitopes.tsv found under $PVAC_OUT/MHC_Class_I"
  fi
  mark_done "pvacseq_${P}"
  info "$P: pVACseq results in $PVAC_OUT"
done
info "done"
