#!/usr/bin/env bash
# =============================================================================
# scripts/lib/common.sh — shared helpers. SOURCE this file; never execute it.
# Callers must have run:  set -Eeuo pipefail
# =============================================================================
[[ -n "${_MELPIPE_COMMON_LOADED:-}" ]] && return 0
_MELPIPE_COMMON_LOADED=1

_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$_LIB_DIR/../.." && pwd)"
export PROJECT_ROOT
# shellcheck disable=SC1091
source "$PROJECT_ROOT/config/config.sh"

: "${FORCE:=0}" "${DRY_RUN:=0}" "${PATIENT_FILTER:=}" "${STEP_NAME:=main}" "${CUR_LOG:=}"
export FORCE DRY_RUN KEEP_INTERMEDIATE PATIENT_FILTER

# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------
_ts()   { date '+%Y-%m-%d %H:%M:%S'; }
_log()  { local lvl="$1"; shift; printf '[%s] [%-5s] [%s] %s\n' "$(_ts)" "$lvl" "$STEP_NAME" "$*" >&2; }
info()  { _log INFO "$@"; }
warn()  { _log WARN "$@"; }
err()   { _log ERROR "$@"; }
die()   { trap - ERR; err "$@"; exit 1; }

_on_err() {
  local ec=$? line="$1" cmd="$2"
  trap - ERR
  err "FAILED (exit ${ec}) at ${BASH_SOURCE[1]:-?}:${line}: ${cmd}"
  if [[ -n "${CUR_LOG:-}" && -s "${CUR_LOG}" ]]; then
    err "---- last lines of ${CUR_LOG} ----"
    tail -n 25 "$CUR_LOG" >&2 || true
  fi
}

_CLEANUP_PATHS=()
register_cleanup() { _CLEANUP_PATHS+=("$@"); }
_on_exit() {
  local ec=$? p
  for p in "${_CLEANUP_PATHS[@]:-}"; do
    if [[ -n "$p" && -e "$p" ]]; then rm -rf -- "$p"; fi
  done
  return "$ec"
}

# init_step NAME — call first thing in every stage script.
init_step() {
  STEP_NAME="$1"
  mkdir -p "$LOG_DIR" "$TMP_DIR" "$MARK_DIR"
  if [[ -z "${MELPIPE_LOGGING:-}" ]]; then      # run standalone: tee to a log file
    local lf; lf="$LOG_DIR/${STEP_NAME}_$(date +%Y%m%d_%H%M%S).log"
    exec > >(tee -a "$lf") 2>&1
    export MELPIPE_LOGGING=1
  fi
  trap '_on_err "$LINENO" "$BASH_COMMAND"' ERR
  trap '_on_exit' EXIT
  load_samplesheet
  info "start (force=$FORCE dry_run=$DRY_RUN patient_filter='${PATIENT_FILTER:-all}')"
}

# ----------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------
trim() { local s="$1"; s="${s#"${s%%[![:space:]]*}"}"; s="${s%"${s##*[![:space:]]}"}"; printf '%s' "$s"; }

require_file() {   # require_file FILE [hint]
  local f="$1" hint="${2:-}"
  if [[ -s "$f" ]]; then return 0; fi
  if (( DRY_RUN )); then warn "[dry-run] missing input (would fail): $f"; return 0; fi
  die "Required file missing or empty: $f ${hint:+— $hint}"
}

require_dir() {
  local d="$1" hint="${2:-}"
  if [[ -d "$d" ]]; then return 0; fi
  if (( DRY_RUN )); then warn "[dry-run] missing dir (would fail): $d"; return 0; fi
  die "Required directory missing: $d ${hint:+— $hint}"
}

free_gb() { mkdir -p "$1"; df -Pk "$1" | awk 'NR==2{printf "%d", $4/1024/1024}'; }
require_free_gb() {   # require_free_gb GB [path]
  local need="$1" path="${2:-$RESULTS_DIR}" have
  have="$(free_gb "$path")"
  if (( have < need )); then
    if (( DRY_RUN )); then warn "[dry-run] only ${have} GB free on $path (need ${need})"; return 0; fi
    die "Need >= ${need} GB free on $path, only ${have} GB available. Run: ./run.sh cleanup"
  fi
  info "disk OK: ${have} GB free (need ${need})"
}

# ----------------------------------------------------------------------------
# Resume machinery.  Each unit of work has a marker; a unit is skipped only if
# its marker AND all of its expected output files exist. --force ignores markers.
# ----------------------------------------------------------------------------
is_done() {   # is_done UNIT [output files...]
  local unit="$1"; shift
  if (( FORCE )); then return 1; fi
  [[ -f "$MARK_DIR/${unit}.done" ]] || return 1
  local f
  for f in "$@"; do [[ -s "$f" ]] || return 1; done
  return 0
}
mark_done() {
  if (( DRY_RUN )); then return 0; fi
  mkdir -p "$MARK_DIR"
  date -Is > "$MARK_DIR/${1}.done"
}
unmark() { rm -f "$MARK_DIR/${1}.done"; }

# Remove intermediates unless KEEP_INTERMEDIATE=1
discard() {
  if (( KEEP_INTERMEDIATE )); then info "keeping intermediates (KEEP_INTERMEDIATE=1): $*"; return 0; fi
  local f
  for f in "$@"; do
    if [[ -e "$f" ]]; then info "cleanup: removing $f"; (( DRY_RUN )) || rm -rf -- "$f"; fi
  done
  return 0
}

# ----------------------------------------------------------------------------
# Samplesheet  (patient,sample,lane,fastq_1,fastq_2,status)
# ----------------------------------------------------------------------------
declare -a SS_PATIENT=() SS_SAMPLE=() SS_LANE=() SS_FQ1=() SS_FQ2=() SS_STATUS=()

load_samplesheet() {
  [[ -f "$SAMPLESHEET" ]] || die "Samplesheet not found: $SAMPLESHEET"
  SS_PATIENT=(); SS_SAMPLE=(); SS_LANE=(); SS_FQ1=(); SS_FQ2=(); SS_STATUS=()
  local line first=1 patient sample lane fq1 fq2 status
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    if [[ -z "${line//[[:space:],]/}" || "$line" =~ ^[[:space:]]*# ]]; then continue; fi
    if (( first )); then
      first=0
      if [[ "${line,,}" == patient,* ]]; then continue; fi
    fi
    IFS=',' read -r patient sample lane fq1 fq2 status <<< "$line"
    patient="$(trim "${patient:-}")"; sample="$(trim "${sample:-}")"; lane="$(trim "${lane:-}")"
    fq1="$(trim "${fq1:-}")"; fq2="$(trim "${fq2:-}")"; status="$(trim "${status:-}")"
    case "${sample,,}" in
      tumor|tumour) sample="Tumor" ;;
      normal)       sample="Normal" ;;
      rna|rnaseq|rna-seq|rna_seq) sample="RNA" ;;
      *) die "Samplesheet: unknown sample type '$sample' (use Tumor, Normal or RNA)" ;;
    esac
    [[ -n "$patient" && -n "$fq1" && -n "$fq2" ]] || die "Samplesheet: incomplete row: $line"
    [[ -z "$lane" ]] && lane=1
    if [[ -n "$PATIENT_FILTER" && "$patient" != "$PATIENT_FILTER" ]]; then continue; fi
    SS_PATIENT+=("$patient"); SS_SAMPLE+=("$sample"); SS_LANE+=("$lane")
    SS_FQ1+=("$fq1"); SS_FQ2+=("$fq2"); SS_STATUS+=("$status")
  done < "$SAMPLESHEET"
  (( ${#SS_PATIENT[@]} > 0 )) || die "Samplesheet has no usable rows (patient filter: '${PATIENT_FILTER:-none}')"
}

patients() { printf '%s\n' "${SS_PATIENT[@]}" | awk '!seen[$0]++'; }

# rows PATIENT SAMPLE -> "lane<TAB>fq1<TAB>fq2" (sorted by lane)
rows() {
  local p="$1" s="$2" i
  for i in "${!SS_PATIENT[@]}"; do
    if [[ "${SS_PATIENT[i]}" == "$p" && "${SS_SAMPLE[i]}" == "$s" ]]; then
      printf '%s\t%s\t%s\n' "${SS_LANE[i]}" "${SS_FQ1[i]}" "${SS_FQ2[i]}"
    fi
  done | sort -k1,1V
  return 0
}
has_sample() { [[ -n "$(rows "$1" "$2")" ]]; }

# ctx PATIENT — set per-patient path variables
ctx() {
  P="$1"
  D_TRIM="$RESULTS_DIR/trimmed/$P"
  D_ALN="$RESULTS_DIR/alignment/$P"
  D_VAR="$RESULTS_DIR/variants/$P"
  D_RNA="$RESULTS_DIR/rna/$P"
  D_HLA="$RESULTS_DIR/hla/$P"
  D_NEO="$RESULTS_DIR/neoantigen/$P"
  D_MRNA="$RESULTS_DIR/mrna/$P"
  D_REP="$RESULTS_DIR/report/$P"
  QC_DIR="$RESULTS_DIR/qc"
  SM_TUMOR="${P}_Tumor"; SM_NORMAL="${P}_Normal"; SM_RNA="${P}_RNA"
  BAM_TUMOR="$D_ALN/${P}_Tumor.md.bam"
  BAM_NORMAL="$D_ALN/${P}_Normal.md.bam"
  BAM_RNA="$D_RNA/star/${P}_RNA.Aligned.sortedByCoord.out.bam"
  QUANT_SF="$D_RNA/salmon/quant.sf"
  VCF_PASS="$D_VAR/${P}.somatic.pass.vcf.gz"
  VCF_NORM="$D_VAR/${P}.somatic.pass.norm.vcf.gz"
  VCF_VEP="$D_VAR/${P}.somatic.pass.norm.vep.vcf.gz"
  HLA_FILE="$D_HLA/${P}.hla_alleles.txt"
  PVAC_INPUT="$D_NEO/${P}.pvac_input.vcf.gz"
  PVAC_OUT="$D_NEO/pvacseq"
  RANKED="$D_NEO/${P}.ranked_epitopes.tsv"
  DESIGN_INPUT="$D_NEO/${P}.design_input.tsv"
}

# downstream_complete PATIENT SAMPLE — true if the consumer of the trimmed FASTQ
# already finished, so fastp/alignment need not be repeated after cleanup.
downstream_complete() {
  local p="$1" s="$2"
  case "$s" in
    Tumor|Normal) is_done "markdup_${p}_${s}" "$RESULTS_DIR/alignment/$p/${p}_${s}.md.bam" ;;
    RNA) is_done "star_${p}" "$RESULTS_DIR/rna/$p/star/${p}_RNA.Aligned.sortedByCoord.out.bam" \
         && is_done "salmon_${p}" "$RESULTS_DIR/rna/$p/salmon/quant.sf" ;;
    *) return 1 ;;
  esac
}

# ----------------------------------------------------------------------------
# Docker wrapper.  Paths are mounted at the SAME location inside the container,
# so host paths can be passed straight to tools.
#   dk IMAGE [-m MEM] [-c CPUS] [-e K=V] [-v H:C[:ro]] [-w DIR] [-i] -- cmd args...
# ----------------------------------------------------------------------------
_MOUNT_ARGS=()
_SEL=()
_covered() {
  local d="$1" m
  for m in "${_SEL[@]:-}"; do
    if [[ -n "$m" && ( "$d" == "$m" || "$d" == "$m"/* ) ]]; then return 0; fi
  done
  return 1
}
build_mounts() {
  _MOUNT_ARGS=(); _SEL=()
  local d i rw=() ro=()
  rw=("$PROJECT_ROOT" "$REF_DIR" "$DATA_DIR" "$RESULTS_DIR" "$LOG_DIR" "$TMP_DIR" "${EXTRA_RW_MOUNTS[@]:-}")
  for i in "${!SS_FQ1[@]}"; do ro+=("$(dirname "${SS_FQ1[i]}")" "$(dirname "${SS_FQ2[i]}")"); done
  ro+=("${EXTRA_RO_MOUNTS[@]:-}")
  if [[ -n "${TARGET_BED:-}" ]]; then ro+=("$(dirname "$TARGET_BED")"); fi
  if [[ -n "${PVAC_LOCAL_IEDB_HOST_DIR:-}" ]]; then ro+=("$PVAC_LOCAL_IEDB_HOST_DIR"); fi
  # shortest paths first so parents win over children
  local sorted
  sorted="$(printf '%s\n' "${rw[@]}" | awk 'NF' | while read -r d; do realpath -m "$d"; done | awk '{print length($0)"\t"$0}' | sort -n | cut -f2- | awk '!s[$0]++')"
  while IFS= read -r d; do
    [[ -n "$d" ]] || continue
    if ! _covered "$d"; then _SEL+=("$d"); _MOUNT_ARGS+=(-v "$d:$d"); fi
  done <<< "$sorted"
  sorted="$(printf '%s\n' "${ro[@]}" | awk 'NF' | while read -r d; do realpath -m "$d"; done | awk '{print length($0)"\t"$0}' | sort -n | cut -f2- | awk '!s[$0]++')"
  while IFS= read -r d; do
    [[ -n "$d" ]] || continue
    if ! _covered "$d"; then _SEL+=("$d"); _MOUNT_ARGS+=(-v "$d:$d:ro"); fi
  done <<< "$sorted"
}

dk() {
  local image="$1"; shift
  local mem="" cpus="" workdir="$PROJECT_ROOT" interactive=0
  local -a extra=()
  while (( $# )); do
    case "$1" in
      --) shift; break ;;
      -m) mem="$2"; shift 2 ;;
      -c) cpus="$2"; shift 2 ;;
      -e) extra+=(-e "$2"); shift 2 ;;
      -v) extra+=(-v "$2"); shift 2 ;;
      -w) workdir="$2"; shift 2 ;;
      -i) interactive=1; shift ;;
      *)  die "dk: unknown option '$1'" ;;
    esac
  done
  if (( ${#_MOUNT_ARGS[@]} == 0 )); then build_mounts; fi
  local -a cmd=(docker run --rm --init
                --user "$(id -u):$(id -g)"
                --label "melpipe_run=${MELPIPE_RUN_ID:-manual}"
                -e HOME=/tmp -w "$workdir" "${_MOUNT_ARGS[@]}")
  if (( interactive )); then cmd+=(-i); fi
  if [[ -n "$mem" ]];  then cmd+=(--memory "$mem"); fi
  if [[ -n "$cpus" ]]; then cmd+=(--cpus "$cpus"); fi
  cmd+=("${extra[@]}" "$image" "$@")
  if (( DRY_RUN )); then
    info "[dry-run] ${cmd[*]}"
    return 0
  fi
  "${cmd[@]}"
}

image_present() { docker image inspect "$1" >/dev/null 2>&1; }

# ensure_image IMAGE — pull (or build, for melpipe/*:local) if absent
ensure_image() {
  local img="$1"
  if (( DRY_RUN )); then return 0; fi
  if image_present "$img"; then return 0; fi
  if [[ -n "${LOCAL_IMAGES[$img]:-}" ]]; then
    die "Local image '$img' not built yet. Run:  ./run.sh setup --components images"
  fi
  info "pulling $img"
  docker pull "$img" || die "Cannot pull $img — edit the tag in config/config.sh (IMG_*)"
}

# gatk_run XMX DOCKER_MEM TOOL [args...]
gatk_run() {
  local xmx="$1" mem="$2"; shift 2
  dk "$IMG_GATK" -m "$mem" -- gatk \
     --java-options "-Xmx${xmx} -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${TMP_DIR}" "$@"
}

# run_logged LOGFILE cmd... — merge stdout+stderr into LOGFILE, show tail on failure
run_logged() {
  local lf="$1"; shift
  mkdir -p "$(dirname "$lf")"
  CUR_LOG="$lf"
  info "log: $lf"
  if (( DRY_RUN )); then "$@"; CUR_LOG=""; return 0; fi
  if ! "$@" >>"$lf" 2>&1; then
    return 1
  fi
  CUR_LOG=""
}

# download URL DEST — resumable, atomic
download() {
  local url="$1" dest="$2"
  if [[ -s "$dest" ]] && (( ! FORCE )); then info "exists: $dest"; return 0; fi
  mkdir -p "$(dirname "$dest")"
  info "download: $url"
  if (( DRY_RUN )); then return 0; fi
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 5 --retry-delay 5 -C - -o "${dest}.partial" "$url"
  else
    wget -c -O "${dest}.partial" "$url"
  fi
  mv -f "${dest}.partial" "$dest"
}

# Helper for FASTQ-derived base names (used by FastQC)
fq_base() { local b; b="$(basename "$1")"; b="${b%.gz}"; b="${b%.fastq}"; b="${b%.fq}"; printf '%s' "$b"; }

# Human-readable duration
fmt_dur() { local s="$1"; printf '%dh%02dm%02ds' $((s/3600)) $((s%3600/60)) $((s%60)); }
