#!/usr/bin/env bash
# Stage alignment (1/2): BWA-MEM (or BWA-MEM2) -> coordinate-sorted BAM per lane -> merged BAM.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "03_bwa_mem2"

ensure_image "$IMG_ALIGN"
case "$ALIGNER" in bwa|bwa-mem2) ;; *) die "ALIGNER must be bwa or bwa-mem2 (got '$ALIGNER')" ;; esac
require_file "$REF_FASTA"
if [[ "$ALIGNER" == "bwa-mem2" ]]; then require_file "${REF_FASTA}.bwt.2bit.64" "prebuilt bwa-mem2 index"; else require_file "${REF_FASTA}.bwt" "run ./run.sh setup"; fi

align_lane() {   # P S lane r1 r2 out_bam
  local P="$1" S="$2" lane="$3" r1="$4" r2="$5" out="$6"
  local rg="@RG\\tID:${P}_${S}_L${lane}\\tSM:${P}_${S}\\tPL:ILLUMINA\\tLB:${P}_${S}_lib1\\tPU:${P}_${S}_L${lane}"
  local part="${out%.bam}.partial.bam" lf="$LOG_DIR/align_${P}_${S}_L${lane}.log"
  CUR_LOG="$lf"
  dk "$IMG_ALIGN" -m "$ALIGN_DOCKER_MEM" -c "$THREADS" \
     -e "RG=$rg" -e "REF=$REF_FASTA" -e "R1=$r1" -e "R2=$r2" -e "OUT=$part" \
     -e "NT=$ALIGN_THREADS" -e "ST=$SORT_THREADS" -e "SM=$SORT_MEM_PER_THREAD" \
     -e "TMPP=$TMP_DIR/sort_${P}_${S}_L${lane}" -e "ALN=$ALIGNER" -- \
     bash -c 'set -euo pipefail
              $ALN mem -t "$NT" -K 100000000 -Y -R "$RG" "$REF" "$R1" "$R2" \
                | samtools sort -@ "$ST" -m "$SM" -T "$TMPP" -O bam -o "$OUT" -
              samtools quickcheck -v "$OUT"' 2>>"$lf"
  CUR_LOG=""
  if (( ! DRY_RUN )); then mv -f "$part" "$out"; fi
}

for P in $(patients); do
  ctx "$P"; mkdir -p "$D_ALN" "$D_TRIM"
  for S in Tumor Normal; do
    if ! has_sample "$P" "$S"; then die "$P: samplesheet has no $S row"; fi
    final="$D_ALN/${P}_${S}.md.bam"; merged="$D_ALN/${P}_${S}.merged.bam"
    if is_done "markdup_${P}_${S}" "$final"; then info "skip (mark-duplicates output exists): $P $S"; continue; fi
    if is_done "merge_${P}_${S}" "$merged"; then info "skip (merged BAM exists): $P $S"; continue; fi
    require_free_gb 50 "$RESULTS_DIR"
    lane_bams=(); trimmed=()
    while IFS=$'\t' read -r lane _ _; do
      t1="$D_TRIM/${P}_${S}_L${lane}_R1.fastq.gz"; t2="$D_TRIM/${P}_${S}_L${lane}_R2.fastq.gz"
      lb="$D_ALN/${P}_${S}_L${lane}.sorted.bam"
      lane_bams+=("$lb"); trimmed+=("$t1" "$t2")
      if is_done "align_${P}_${S}_L${lane}" "$lb"; then info "skip (done): lane $lane"; continue; fi
      require_file "$t1" "run stage qc first (02_fastp)"; require_file "$t2"
      info "aligning $P $S lane $lane with $ALIGNER"
      align_lane "$P" "$S" "$lane" "$t1" "$t2" "$lb"
      mark_done "align_${P}_${S}_L${lane}"
    done < <(rows "$P" "$S")

    if (( ${#lane_bams[@]} == 1 )); then
      info "single lane: using lane BAM as merged BAM"
      if (( ! DRY_RUN )); then mv -f "${lane_bams[0]}" "$merged"; fi
    else
      info "merging ${#lane_bams[@]} lane BAMs"
      run_logged "$LOG_DIR/merge_${P}_${S}.log" \
        dk "$IMG_ALIGN" -m 8g -c 8 -- samtools merge -@ 8 -f "${merged%.bam}.partial.bam" "${lane_bams[@]}"
      if (( ! DRY_RUN )); then mv -f "${merged%.bam}.partial.bam" "$merged"; fi
      discard "${lane_bams[@]}"
    fi
    mark_done "merge_${P}_${S}"
    if (( CLEAN_TRIMMED_FASTQ )); then discard "${trimmed[@]}"; fi
  done
done
info "done"
