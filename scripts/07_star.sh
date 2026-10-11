#!/usr/bin/env bash
# Stage RNA (1/2): STAR alignment of tumor RNA-seq (sparse-SA index to fit in 28 GB RAM).
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "07_star"

ensure_image "$IMG_STAR"
require_file "$STAR_INDEX_DIR/SA" "run ./run.sh setup --components star"

for P in $(patients); do
  ctx "$P"
  has_sample "$P" RNA || { warn "$P: no RNA row, skipping"; continue; }
  S="$D_RNA/star"; mkdir -p "$S"
  if is_done "star_${P}" "$BAM_RNA" "$BAM_RNA.bai"; then info "skip (done): $P"; else
    r1=(); r2=()
    while IFS=$'\t' read -r lane _ _; do
      a="$D_TRIM/${P}_RNA_L${lane}_R1.fastq.gz"; b="$D_TRIM/${P}_RNA_L${lane}_R2.fastq.gz"
      require_file "$a" "run stage qc first"; require_file "$b"
      r1+=("$a"); r2+=("$b")
    done < <(rows "$P" RNA)
    require_free_gb 30 "$RESULTS_DIR"
    in1="$(IFS=,; echo "${r1[*]}")"; in2="$(IFS=,; echo "${r2[*]}")"
    tmpd="$TMP_DIR/star_${P}"; rm -rf "$tmpd"; register_cleanup "$tmpd"
    info "STAR $P (index ~18 GB in RAM; ~1 h for typical RNA-seq)"
    run_logged "$LOG_DIR/star_${P}.log" \
      dk "$IMG_STAR" -m "$STAR_DOCKER_MEM" -c "$STAR_THREADS" -- STAR \
        --runThreadN "$STAR_THREADS" --genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
        --readFilesIn "$in1" "$in2" --readFilesCommand zcat \
        --outFileNamePrefix "$S/${P}_RNA." --outTmpDir "$tmpd" \
        --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 --limitBAMsortRAM "$STAR_SORT_RAM" \
        --outSAMattrRGline "ID:${P}_RNA" "SM:${P}_RNA" PL:ILLUMINA \
        --outSAMattributes NH HI AS nM NM MD jM jI MC ch --outSAMstrandField intronMotif \
        --outFilterMultimapNmax 10 --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 --twopassMode "$STAR_TWOPASS" --quantMode GeneCounts
    run_logged "$LOG_DIR/star_index_bam_${P}.log" \
      dk "$IMG_STAR" -m 4g -c 8 -- samtools index -@ 8 "$BAM_RNA"
    mark_done "star_${P}"
  fi
  # free trimmed RNA reads once both STAR and salmon are finished
  if (( CLEAN_TRIMMED_FASTQ )) && is_done "star_${P}" "$BAM_RNA" && is_done "salmon_${P}" "$QUANT_SF"; then
    discard "$D_TRIM"/${P}_RNA_L*_R[12].fastq.gz
  fi
done
info "done"
