#!/usr/bin/env bash
# Stage alignment (2/2): GATK MarkDuplicates, index, flagstat/idxstats, mosdepth coverage.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "04_markdup"

ensure_image "$IMG_GATK"; ensure_image "$IMG_ALIGN"
MOSD="$RESULTS_DIR/qc/mosdepth"; mkdir -p "$MOSD"

for P in $(patients); do
  ctx "$P"; mkdir -p "$D_ALN"
  for S in Tumor Normal; do
    merged="$D_ALN/${P}_${S}.merged.bam"; final="$D_ALN/${P}_${S}.md.bam"
    if is_done "markdup_${P}_${S}" "$final" "$final.bai"; then info "skip (done): $P $S"; continue; fi
    require_file "$merged" "run stage alignment step 03 first"
    require_free_gb 25 "$RESULTS_DIR"
    info "MarkDuplicates $P $S"
    part="${final%.bam}.partial.bam"
    run_logged "$LOG_DIR/markdup_${P}_${S}.log" \
      gatk_run "$MARKDUP_XMX" "$MARKDUP_DOCKER_MEM" MarkDuplicates \
        -I "$merged" -O "$part" -M "$D_ALN/${P}_${S}.markdup_metrics.txt" \
        --ASSUME_SORT_ORDER coordinate --CREATE_INDEX false --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE "$OPTICAL_DUP_DISTANCE" --COMPRESSION_LEVEL 5 --TMP_DIR "$TMP_DIR"
    if (( ! DRY_RUN )); then mv -f "$part" "$final"; fi
    run_logged "$LOG_DIR/bamstats_${P}_${S}.log" \
      dk "$IMG_ALIGN" -m 6g -c 8 -e "BAM=$final" -e "FS=$D_ALN/${P}_${S}.flagstat.txt" -e "IS=$D_ALN/${P}_${S}.idxstats.txt" -- \
        bash -c 'set -euo pipefail; samtools index -@ 8 "$BAM"; samtools flagstat -@ 8 "$BAM" > "$FS"; samtools idxstats "$BAM" > "$IS"'

    # coverage: quantized COVERED/LOW bed (also used to derive exome targets) + summary
    byargs=(); [[ -n "$TARGET_BED" ]] && byargs=(--by "$TARGET_BED")
    run_logged "$LOG_DIR/mosdepth_${P}_${S}.log" \
      dk "$IMG_ALIGN" -m 6g -c "$MOSDEPTH_THREADS" -e MOSDEPTH_Q0=LOW -e MOSDEPTH_Q1=COVERED -- \
        mosdepth -t "$MOSDEPTH_THREADS" -n --quantize "0:${MIN_COVERED_DEPTH}:" "${byargs[@]}" "$MOSD/${P}_${S}" "$final"
    discard "$merged"
    mark_done "markdup_${P}_${S}"
  done
done
info "done"
