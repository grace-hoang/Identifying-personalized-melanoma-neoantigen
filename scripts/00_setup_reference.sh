#!/usr/bin/env bash
# Download and index every reference resource. Fully resumable.
#   SETUP_COMPONENTS="genome gatk gencode bwa star salmon vep"   (default: all)
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "00_setup_reference"

COMPONENTS="${SETUP_COMPONENTS:-all}"
want() { [[ "$COMPONENTS" == "all" || " $COMPONENTS " == *" $1 "* ]]; }

mkdir -p "$REF_DIR"/{genome,gatk,annotation,index,vep} "$TMP_DIR"
require_free_gb 20 "$REF_DIR"

# ---------------------------------------------------------------- genome ----
if want genome; then
  download "$URL_BROAD_HG38/Homo_sapiens_assembly38.fasta"     "$REF_FASTA"
  download "$URL_BROAD_HG38/Homo_sapiens_assembly38.fasta.fai" "$REF_FAI"
  download "$URL_BROAD_HG38/Homo_sapiens_assembly38.dict"      "$REF_DICT"
  if [[ ! -s "$PRIMARY_CONTIGS_LIST" ]] || (( FORCE )); then
    awk -v re="$PRIMARY_CONTIG_REGEX" '$1 ~ re {print $1}' "$REF_FAI" > "$PRIMARY_CONTIGS_LIST"
    info "primary contigs: $(wc -l < "$PRIMARY_CONTIGS_LIST")"
  fi
fi

# ------------------------------------------------------------------ gatk ----
if want gatk; then
  for f in af-only-gnomad.hg38.vcf.gz 1000g_pon.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz; do
    download "$URL_GATK_SOMATIC/$f"     "$REF_DIR/gatk/$f"
    download "$URL_GATK_SOMATIC/$f.tbi" "$REF_DIR/gatk/$f.tbi"
  done
fi

# --------------------------------------------------------------- gencode ----
if want gencode || want star || want salmon; then
  download "$URL_GENCODE/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" "$GENCODE_GTF.gz"
  download "$URL_GENCODE/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz" "$GENCODE_TX_FASTA"
  download "$URL_GENCODE/GRCh38.primary_assembly.genome.fa.gz"          "$PRIMARY_FASTA.gz"
  if [[ ! -s "$GENCODE_GTF" ]]; then
    info "decompressing GTF"
    if (( ! DRY_RUN )); then
      gzip -dc "$GENCODE_GTF.gz" > "$GENCODE_GTF.partial"
      mv -f "$GENCODE_GTF.partial" "$GENCODE_GTF"
    fi
  fi
  if [[ ! -s "$TX2GENE" ]]; then
    info "building tx2gene table"
    if (( ! DRY_RUN )); then
      awk -F'\t' 'BEGIN{OFS="\t"; print "transcript_id","gene_id","gene_name"}
        $3=="transcript"{tid="";gid="";gn="";n=split($9,a,";")
          for(i=1;i<=n;i++){gsub(/^ +| +$/,"",a[i])
            if(a[i]~/^transcript_id/){tid=a[i];gsub(/transcript_id "|"/,"",tid)}
            else if(a[i]~/^gene_id/){gid=a[i];gsub(/gene_id "|"/,"",gid)}
            else if(a[i]~/^gene_name/){gn=a[i];gsub(/gene_name "|"/,"",gn)}}
          print tid,gid,gn}' "$GENCODE_GTF" > "$TX2GENE.partial"
      mv -f "$TX2GENE.partial" "$TX2GENE"
    fi
  fi
fi

# ------------------------------------------------------------------- bwa ----
if want bwa; then
  ensure_image "$IMG_ALIGN"
  if [[ "$ALIGNER" == "bwa-mem2" ]]; then
    idx="${REF_FASTA}.bwt.2bit.64"
    mem_gb=$(awk '/MemTotal/{printf "%d",$2/1024/1024}' /proc/meminfo)
    if [[ ! -s "$idx" ]] && (( mem_gb < 80 )); then
      die "bwa-mem2 index construction needs ~28 GB RAM per Gbp (~87 GB for hg38); this machine has ${mem_gb} GB. Use ALIGNER=bwa, or copy a prebuilt bwa-mem2 index next to $REF_FASTA."
    fi
  fi
  if [[ ! -s "${REF_FASTA}.bwt" ]] || (( FORCE )); then
    info "bwa index (single-threaded, ~1 h, ~6 GB RAM)"
    CUR_LOG="$LOG_DIR/bwa_index.log"
    dk "$IMG_ALIGN" -m 12g -- bwa index -a bwtsw "$REF_FASTA" >>"$CUR_LOG" 2>&1
    CUR_LOG=""
  else info "bwa index present"; fi
  if [[ "$ALIGNER" == "bwa-mem2" && ! -s "${REF_FASTA}.bwt.2bit.64" ]]; then
    CUR_LOG="$LOG_DIR/bwa_mem2_index.log"
    dk "$IMG_ALIGN" -- bwa-mem2 index "$REF_FASTA" >>"$CUR_LOG" 2>&1
    CUR_LOG=""
  fi
fi

# ------------------------------------------------------------------ star ----
if want star; then
  ensure_image "$IMG_STAR"
  if [[ ! -s "$PRIMARY_FASTA" ]]; then
    info "decompressing primary assembly FASTA"
    (( DRY_RUN )) || { gzip -dc "$PRIMARY_FASTA.gz" > "$PRIMARY_FASTA.partial"; mv -f "$PRIMARY_FASTA.partial" "$PRIMARY_FASTA"; }
  fi
  if [[ ! -s "$STAR_INDEX_DIR/SA" ]] || (( FORCE )); then
    require_free_gb 40 "$REF_DIR"
    info "STAR genomeGenerate (sparse SA D=${STAR_SA_SPARSE}); needs ~22-26 GB RAM, 30-60 min"
    mkdir -p "$STAR_INDEX_DIR"
    rm -rf "$TMP_DIR/star_genome_tmp"
    CUR_LOG="$LOG_DIR/star_index.log"
    dk "$IMG_STAR" -m 26g -c "$STAR_THREADS" -- STAR --runMode genomeGenerate \
       --runThreadN "$STAR_THREADS" --genomeDir "$STAR_INDEX_DIR" \
       --genomeFastaFiles "$PRIMARY_FASTA" --sjdbGTFfile "$GENCODE_GTF" \
       --sjdbOverhang "$STAR_SJDB_OVERHANG" --genomeSAsparseD "$STAR_SA_SPARSE" \
       --genomeSAindexNbases 14 --limitGenomeGenerateRAM 25000000000 \
       --outTmpDir "$TMP_DIR/star_genome_tmp" --outFileNamePrefix "$STAR_INDEX_DIR/" >>"$CUR_LOG" 2>&1
    CUR_LOG=""
  else info "STAR index present"; fi
fi

# ---------------------------------------------------------------- salmon ----
if want salmon; then
  ensure_image "$IMG_SALMON"
  if [[ ! -s "$SALMON_INDEX_DIR/versionInfo.json" ]] || (( FORCE )); then
    require_free_gb 25 "$REF_DIR"
    info "salmon index with genome decoys (selective alignment)"
    gent="$REF_DIR/index/gentrome.gencode${GENCODE_RELEASE}.fa.gz"
    dec="$REF_DIR/index/decoys.gencode${GENCODE_RELEASE}.txt"
    if (( ! DRY_RUN )); then
      gzip -dc "$PRIMARY_FASTA.gz" | grep '^>' | cut -d' ' -f1 | sed 's/^>//' > "$dec"
      cat "$GENCODE_TX_FASTA" "$PRIMARY_FASTA.gz" > "$gent"
    fi
    rm -rf "$SALMON_INDEX_DIR"
    CUR_LOG="$LOG_DIR/salmon_index.log"
    dk "$IMG_SALMON" -m "$SALMON_DOCKER_MEM" -c "$SALMON_THREADS" -- salmon index \
       -t "$gent" -d "$dec" -i "$SALMON_INDEX_DIR" -p "$SALMON_THREADS" -k 31 --gencode >>"$CUR_LOG" 2>&1
    CUR_LOG=""
    discard "$gent"
  else info "salmon index present"; fi
fi

# ------------------------------------------------------------------- vep ----
if want vep; then
  ensure_image "$IMG_PVAC"
  mkdir -p "$VEP_CACHE_DIR" "$VEP_PLUGIN_DIR"
  if [[ ! -d "$VEP_CACHE_DIR/homo_sapiens/${ENSEMBL_RELEASE}_GRCh38" ]] || (( FORCE )); then
    require_free_gb 45 "$REF_DIR"
    tarball="$REF_DIR/vep/$(basename "$URL_VEP_CACHE")"
    download "$URL_VEP_CACHE" "$tarball"
    info "extracting VEP cache (~15-20 GB)"
    (( DRY_RUN )) || tar -xzf "$tarball" -C "$VEP_CACHE_DIR"
    discard "$tarball"
  else info "VEP cache present"; fi
  if [[ ! -s "$VEP_PLUGIN_DIR/Wildtype.pm" || ! -s "$VEP_PLUGIN_DIR/Frameshift.pm" ]] || (( FORCE )); then
    info "installing pVACseq VEP plugins (Frameshift, Wildtype)"
    dk "$IMG_PVAC" -- pvacseq install_vep_plugin "$VEP_PLUGIN_DIR"
  else info "VEP plugins present"; fi
fi

info "reference setup finished. Disk used by reference/: $(du -sh "$REF_DIR" 2>/dev/null | cut -f1)"
