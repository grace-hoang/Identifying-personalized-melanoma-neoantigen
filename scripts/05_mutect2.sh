#!/usr/bin/env bash
# Stage variants (1/2): somatic SNV/indel calling with GATK Mutect2 (tumor vs matched normal),
# scattered over exome targets, then contamination estimation + FilterMutectCalls.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "05_mutect2"

ensure_image "$IMG_GATK"; ensure_image "$IMG_ALIGN"
for f in "$REF_FASTA" "$REF_FAI" "$REF_DICT" "$GERMLINE_RESOURCE" "$PON" "$CONTAM_VARIANTS"; do require_file "$f" "run ./run.sh setup"; done

for P in $(patients); do
  ctx "$P"; mkdir -p "$D_VAR"
  if is_done "mutect2_${P}" "$VCF_PASS" "$VCF_PASS.tbi"; then info "skip (done): $P"; continue; fi
  require_file "$BAM_TUMOR" "run alignment stage"; require_file "$BAM_NORMAL" "run alignment stage"
  require_free_gb 20 "$RESULTS_DIR"
  W="$D_VAR/mutect2_work"; mkdir -p "$W/scatter"
  TARGETS="$D_VAR/${P}.targets.bed"

  # ---- 1. target intervals ---------------------------------------------------
  if [[ -n "$TARGET_BED" ]]; then src="$TARGET_BED"; flt="all"
  else src="$RESULTS_DIR/qc/mosdepth/${P}_Normal.quantized.bed.gz"; flt="covered"; fi
  if is_done "targets_${P}" "$TARGETS"; then info "skip (done): targets"; else
    require_file "$src" "needs 04_markdup mosdepth output or TARGET_BED"
    info "building calling intervals from $(basename "$src") (filter: $flt)"
    run_logged "$LOG_DIR/targets_${P}.log" \
      dk "$IMG_ALIGN" -m 4g -e "SRC=$src" -e "FLT=$flt" -e "RE=$PRIMARY_CONTIG_REGEX" -e "FAI=$REF_FAI" \
         -e "PAD=$TARGET_PAD" -e "OUT=$TARGETS.partial" -- bash -c '
        set -euo pipefail
        { if [[ "$SRC" == *.gz ]]; then zcat "$SRC"; else cat "$SRC"; fi; } \
          | awk -F"\t" -v re="$RE" -v flt="$FLT" '"'"'$1 ~ re && (flt=="all" || $4=="COVERED") {print $1"\t"$2"\t"$3}'"'"' \
          | sort -k1,1 -k2,2n | bedtools merge -d 150 -i - | awk '"'"'($3-$2)>=20'"'"' \
          | bedtools slop -i - -g "$FAI" -b "$PAD" | sort -k1,1 -k2,2n | bedtools merge -i - > "$OUT"'
    if (( ! DRY_RUN )); then
      mv -f "$TARGETS.partial" "$TARGETS"
      info "target intervals: $(wc -l < "$TARGETS") regions, $(awk '{s+=$3-$2} END{printf "%.1f", s/1e6}' "$TARGETS") Mb"
    fi
    mark_done "targets_${P}"
  fi

  # ---- 2. scatter ------------------------------------------------------------
  if ! is_done "scatter_${P}" "$W/scatter/0000-scattered.interval_list"; then
    rm -f "$W"/scatter/*.interval_list
    run_logged "$LOG_DIR/scatter_${P}.log" \
      gatk_run "$GATK_POST_XMX" "$GATK_POST_DOCKER_MEM" SplitIntervals -R "$REF_FASTA" -L "$TARGETS" \
        --scatter-count "$MUTECT_SCATTER" -O "$W/scatter"
    mark_done "scatter_${P}"
  fi

  # ---- 3. Mutect2 shards in parallel ----------------------------------------
  run_shard() {   # idx interval_list
    local idx="$1" il="$2"
    local vcf="$W/shard_${idx}.vcf.gz" f1r2="$W/f1r2_${idx}.tar.gz" unit="mutect2_shard_${P}_${idx}"
    if is_done "$unit" "$vcf" "$vcf.stats" "$f1r2"; then info "skip (done): shard $idx"; return 0; fi
    info "Mutect2 shard $idx"
    gatk_run "$MUTECT_XMX" "$MUTECT_DOCKER_MEM" Mutect2 -R "$REF_FASTA" \
      -I "$BAM_TUMOR" -I "$BAM_NORMAL" -normal "$SM_NORMAL" \
      --germline-resource "$GERMLINE_RESOURCE" --panel-of-normals "$PON" \
      --f1r2-tar-gz "$f1r2" -L "$il" --native-pair-hmm-threads "$PAIRHMM_THREADS" \
      --tmp-dir "$TMP_DIR" -O "$vcf" >>"$LOG_DIR/mutect2_${P}_shard${idx}.log" 2>&1
    mark_done "$unit"
  }
  shards=("$W"/scatter/*.interval_list)
  if (( DRY_RUN )); then shards=("$W/scatter/0000-scattered.interval_list"); fi
  if (( ! DRY_RUN )) && [[ ! -e "${shards[0]}" ]]; then die "no scatter interval lists in $W/scatter"; fi
  pids=(); n=0
  for il in "${shards[@]}"; do
    idx="$(printf '%04d' "$n")"; n=$((n + 1))
    ( run_shard "$idx" "$il" ) & pids+=($!)
  done
  fail=0
  for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
  if (( fail )); then die "one or more Mutect2 shards failed — see $LOG_DIR/mutect2_${P}_shard*.log (re-run to resume)"; fi

  # ---- 4. merge, contamination, filter --------------------------------------
  vin=(); sin=(); fin=()
  for ((k = 0; k < n; k++)); do
    idx="$(printf '%04d' "$k")"
    vin+=(-I "$W/shard_${idx}.vcf.gz"); sin+=(-stats "$W/shard_${idx}.vcf.gz.stats"); fin+=(-I "$W/f1r2_${idx}.tar.gz")
  done
  UNF="$D_VAR/${P}.somatic.unfiltered.vcf.gz"; FILT="$D_VAR/${P}.somatic.filtered.vcf.gz"
  L="$LOG_DIR/mutect2_post_${P}.log"
  post() { gatk_run "$GATK_POST_XMX" "$GATK_POST_DOCKER_MEM" "$@"; }
  run_logged "$L" post MergeVcfs "${vin[@]}" -O "$UNF"
  run_logged "$L" post MergeMutectStats "${sin[@]}" -O "$UNF.stats"
  run_logged "$L" post LearnReadOrientationModel "${fin[@]}" -O "$D_VAR/${P}.read-orientation-model.tar.gz"
  for T in Tumor Normal; do
    bam="$D_ALN/${P}_${T}.md.bam"
    run_logged "$L" post GetPileupSummaries -I "$bam" -V "$CONTAM_VARIANTS" -L "$CONTAM_VARIANTS" -L "$TARGETS" \
      --interval-set-rule INTERSECTION -R "$REF_FASTA" -O "$D_VAR/${P}.${T}.pileups.table"
  done
  run_logged "$L" post CalculateContamination -I "$D_VAR/${P}.Tumor.pileups.table" -matched "$D_VAR/${P}.Normal.pileups.table" \
    -O "$D_VAR/${P}.contamination.table" --tumor-segmentation "$D_VAR/${P}.segments.table"
  run_logged "$L" post FilterMutectCalls -R "$REF_FASTA" -V "$UNF" --stats "$UNF.stats" \
    --contamination-table "$D_VAR/${P}.contamination.table" --tumor-segmentation "$D_VAR/${P}.segments.table" \
    --ob-priors "$D_VAR/${P}.read-orientation-model.tar.gz" -O "$FILT"
  run_logged "$L" post SelectVariants -R "$REF_FASTA" -V "$FILT" --exclude-filtered -O "$VCF_PASS"
  discard "$W"
  mark_done "mutect2_${P}"
  info "$P: PASS variants -> $VCF_PASS"
done
info "done"
