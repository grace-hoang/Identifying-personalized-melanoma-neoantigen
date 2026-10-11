#!/usr/bin/env bash
# Stage HLA: class-I typing (HLA-A/B/C) from the matched NORMAL exome with OptiType.
# MHC-region + unmapped reads (with their mates) are extracted from the BAM first to keep it fast.
# A manual override (config/hla_override.tsv) skips typing for listed patients.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "09_hla"

ensure_image "$IMG_ALIGN"; ensure_image "$IMG_OPTITYPE"; ensure_image "$IMG_PY"

for P in $(patients); do
  ctx "$P"; mkdir -p "$D_HLA"
  if is_done "hla_${P}" "$HLA_FILE"; then info "skip (done): $P ($(cat "$HLA_FILE" 2>/dev/null))"; continue; fi

  # ---- manual override ----------------------------------------------------
  if [[ -s "$HLA_OVERRIDE_FILE" ]] && grep -qE "^${P}[[:space:]]" "$HLA_OVERRIDE_FILE"; then
    info "$P: using manual HLA override from $HLA_OVERRIDE_FILE"
    run_logged "$LOG_DIR/hla_parse_${P}.log" \
      dk "$IMG_PY" -m 1g -- python3 "$PROJECT_ROOT/scripts/py/hla_parse.py" \
        --patient "$P" --override "$HLA_OVERRIDE_FILE" --out-txt "$HLA_FILE" --out-tsv "$D_HLA/${P}.hla_alleles.tsv"
    mark_done "hla_${P}"; continue
  fi

  # ---- extract MHC reads from normal BAM -----------------------------------
  require_file "$BAM_NORMAL" "run alignment stage"
  require_free_gb 10 "$RESULTS_DIR"
  R1="$D_HLA/${P}.mhc_R1.fastq.gz"; R2="$D_HLA/${P}.mhc_R2.fastq.gz"
  if ! is_done "hla_extract_${P}" "$R1" "$R2"; then
    info "$P: extracting MHC-region reads"
    run_logged "$LOG_DIR/hla_extract_${P}.log" \
      dk "$IMG_ALIGN" -m 12g -c 8 -e "BAM=$BAM_NORMAL" -e "FAI=$REF_FAI" -e "REGION=$HLA_MHC_REGION" \
         -e "OUTDIR=$D_HLA" -e "P=$P" -e "TMPD=$TMP_DIR/hla_${P}" -- bash -c '
        set -euo pipefail
        mkdir -p "$TMPD"
        chrom="${REGION%%:*}"; rng="${REGION#*:}"; s="${rng%-*}"; e="${rng#*-}"
        { printf "%s\t%s\t%s\n" "$chrom" "$((s-1))" "$e"
          awk -F"\t" '"'"'$1 ~ /^HLA-/ {print $1"\t0\t"$2}'"'"' "$FAI"; } > "$TMPD/mhc.bed"
        samtools view -@ 8 -b -L "$TMPD/mhc.bed" "$BAM" > "$TMPD/in_region.bam"
        samtools view -@ 8 -b -f 4 "$BAM" > "$TMPD/unmapped.bam"
        samtools merge -@ 8 -f "$TMPD/seed.bam" "$TMPD/in_region.bam" "$TMPD/unmapped.bam"
        samtools view "$TMPD/seed.bam" | cut -f1 | sort -u > "$TMPD/names.txt"
        samtools view -@ 8 -b -N "$TMPD/names.txt" "$BAM" \
          | samtools sort -n -@ 4 -m 1G -T "$TMPD/s" -O bam - \
          | samtools fastq -@ 4 -n -1 "$OUTDIR/${P}.mhc_R1.fastq.gz" -2 "$OUTDIR/${P}.mhc_R2.fastq.gz" -0 /dev/null -s /dev/null -
        rm -rf "$TMPD"'
    mark_done "hla_extract_${P}"
  fi

  # ---- OptiType --------------------------------------------------------------
  OT="$D_HLA/optitype"; mkdir -p "$OT"
  cfg="$D_HLA/optitype_config.ini"
  if (( ! DRY_RUN )); then
    cat > "$cfg" <<CFGEOF
[mapping]
razers3=razers3
threads=${HLA_OPTITYPE_THREADS}

[ilp]
solver=glpk
threads=1

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
CFGEOF
  fi
  info "$P: OptiType (DNA mode)"
  run_logged "$LOG_DIR/optitype_${P}.log" \
    dk "$IMG_OPTITYPE" -m 12g -c "$HLA_OPTITYPE_THREADS" -- \
      OptiTypePipeline.py -i "$R1" "$R2" --dna -v -c "$cfg" -o "$OT" -p "$P"
  res="$(ls -1 "$OT"/*/"${P}"_result.tsv 2>/dev/null | head -n 1 || true)"
  if [[ -z "$res" ]] && (( ! DRY_RUN )); then die "OptiType produced no result table in $OT"; fi
  run_logged "$LOG_DIR/hla_parse_${P}.log" \
    dk "$IMG_PY" -m 1g -- python3 "$PROJECT_ROOT/scripts/py/hla_parse.py" \
      --patient "$P" --optitype "${res:-/nonexistent}" --out-txt "$HLA_FILE" --out-tsv "$D_HLA/${P}.hla_alleles.tsv"
  info "$P: HLA class I = $(cat "$HLA_FILE" 2>/dev/null || echo '?')"
  discard "$R1" "$R2"
  mark_done "hla_${P}"
done
info "done"
