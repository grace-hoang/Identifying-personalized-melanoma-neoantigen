#!/usr/bin/env bash
# Pre-flight checks: system, Docker, disk, samplesheet, FASTQ, images, reference.
#   PREFLIGHT_STAGES="qc align variants rna hla neoantigen report"  (default all)
#   PREFLIGHT_DEEP=1  also run gzip -t on every FASTQ (slow)
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "00_preflight"

STAGES="${PREFLIGHT_STAGES:-qc align variants rna hla neoantigen report}"
wants() { [[ " $STAGES " == *" $1 "* ]]; }
ERRORS=0; WARNINGS=0
ok()   { info "OK    $*"; }
note() { warn "WARN  $*"; WARNINGS=$((WARNINGS + 1)); }
bad()  { err  "FAIL  $*"; ERRORS=$((ERRORS + 1)); }

# ---- system ---------------------------------------------------------------
if [[ "${BASH_VERSINFO[0]}" -lt 4 ]]; then bad "bash >= 4 required"; else ok "bash ${BASH_VERSION}"; fi
cores=$(nproc); mem_gb=$(awk '/MemTotal/{printf "%d",$2/1024/1024}' /proc/meminfo)
swap_gb=$(awk '/SwapTotal/{printf "%d",$2/1024/1024}' /proc/meminfo)
ok "cores=$cores  RAM=${mem_gb} GB  swap=${swap_gb} GB"
(( cores >= THREADS )) || note "only $cores cores but THREADS=$THREADS"
if (( mem_gb < 24 )); then note "RAM ${mem_gb} GB < 24 GB: STAR/Mutect2 stages may be OOM-killed. Raise WSL2 memory in %UserProfile%\\.wslconfig"; fi
if grep -qi microsoft /proc/version 2>/dev/null; then
  ok "WSL2 detected"
  case "$PROJECT_ROOT" in /mnt/*) note "project is on a Windows drive ($PROJECT_ROOT): I/O will be very slow. Move it under /home." ;; esac
fi
for c in awk sed gzip tar curl realpath flock; do
  command -v "$c" >/dev/null 2>&1 && ok "found $c" || note "missing host tool: $c"
done

# ---- docker ---------------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then bad "docker not found in PATH"
elif ! docker info >/dev/null 2>&1; then bad "docker daemon not reachable (is Docker running / user in 'docker' group?)"
else ok "docker $(docker version --format '{{.Server.Version}}' 2>/dev/null || echo '?')"; fi

# ---- disk -----------------------------------------------------------------
have=$(free_gb "$PROJECT_ROOT")
if (( have < 30 )); then bad "only ${have} GB free on project disk"
elif (( have < MIN_FREE_GB_START )); then note "${have} GB free (< ${MIN_FREE_GB_START} GB recommended)"
else ok "${have} GB free on project disk"; fi

# ---- samplesheet ----------------------------------------------------------
ok "samplesheet: $SAMPLESHEET (${#SS_PATIENT[@]} rows)"
for p in $(patients); do
  has_sample "$p" Tumor  || { if wants align || wants variants; then bad "$p: no Tumor row";  fi; }
  has_sample "$p" Normal || { if wants align || wants variants; then bad "$p: no Normal row"; fi; }
  has_sample "$p" RNA    || { if wants rna; then bad "$p: no RNA row (needed for rna/neoantigen stages)"; fi; }
done
first_id() { ( set +o pipefail; gzip -dc -- "$1" 2>/dev/null | head -n 1 ) || true; }
norm_id()  { local x="${1%% *}"; x="${x%/[12]}"; printf '%s' "$x"; }
for i in "${!SS_PATIENT[@]}"; do
  tag="${SS_PATIENT[i]}/${SS_SAMPLE[i]}/L${SS_LANE[i]}"
  f1="${SS_FQ1[i]}"; f2="${SS_FQ2[i]}"
  if [[ ! -r "$f1" || ! -r "$f2" ]]; then bad "$tag: FASTQ not readable ($f1 / $f2)"; continue; fi
  h1="$(first_id "$f1")"; h2="$(first_id "$f2")"
  if [[ "${h1:0:1}" != "@" || "${h2:0:1}" != "@" ]]; then bad "$tag: not a valid gzip FASTQ"; continue; fi
  if [[ "$(norm_id "$h1")" != "$(norm_id "$h2")" ]]; then note "$tag: first read names of R1/R2 differ ('$h1' vs '$h2')"; fi
  if [[ "${PREFLIGHT_DEEP:-0}" == 1 ]]; then
    gzip -t "$f1" && gzip -t "$f2" && ok "$tag: gzip integrity" || bad "$tag: corrupt gzip"
  fi
  ok "$tag: FASTQ pair readable ($(du -h "$f1" | cut -f1) + $(du -h "$f2" | cut -f1))"
done

# ---- images ---------------------------------------------------------------
declare -a need_imgs=()
wants qc         && need_imgs+=("$IMG_FASTQC" "$IMG_FASTP")
wants align      && need_imgs+=("$IMG_ALIGN" "$IMG_GATK")
wants variants   && need_imgs+=("$IMG_ALIGN" "$IMG_GATK" "$IMG_BCFTOOLS" "$IMG_VEP")
wants rna        && need_imgs+=("$IMG_STAR" "$IMG_SALMON" "$IMG_PY")
wants hla        && need_imgs+=("$IMG_OPTITYPE" "$IMG_ALIGN" "$IMG_PY")
wants neoantigen && need_imgs+=("$IMG_PVAC" "$IMG_BAMRC" "$IMG_BCFTOOLS" "$IMG_PY")
wants report     && need_imgs+=("$IMG_PY" "$IMG_FASTQC")
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
  for img in $(printf '%s\n' "${need_imgs[@]:-}" | awk 'NF && !s[$0]++'); do
    if image_present "$img"; then ok "image $img"; else bad "image missing: $img  -> ./run.sh setup --components images"; fi
  done
fi

# ---- reference ------------------------------------------------------------
chk() { local f="$1" what="$2"; if [[ -s "$f" ]]; then ok "ref: $what"; else bad "ref missing: $what ($f) -> ./run.sh setup"; fi; }
if wants align || wants variants; then
  chk "$REF_FASTA" "genome FASTA"; chk "$REF_FAI" "FASTA .fai"; chk "$REF_DICT" "sequence dictionary"
  if [[ "$ALIGNER" == "bwa-mem2" ]]; then chk "${REF_FASTA}.bwt.2bit.64" "bwa-mem2 index"; else chk "${REF_FASTA}.bwt" "bwa index"; fi
fi
if wants variants; then
  chk "$PRIMARY_CONTIGS_LIST" "primary contig list"
  chk "$GERMLINE_RESOURCE" "gnomAD af-only"; chk "${GERMLINE_RESOURCE}.tbi" "gnomAD index"
  chk "$PON" "panel of normals"; chk "${PON}.tbi" "PoN index"
  chk "$CONTAM_VARIANTS" "common variants (contamination)"; chk "${CONTAM_VARIANTS}.tbi" "common variants index"
  [[ -d "$VEP_CACHE_DIR/homo_sapiens/${ENSEMBL_RELEASE}_GRCh38" ]] && ok "ref: VEP cache" || bad "ref missing: VEP cache -> ./run.sh setup"
  chk "$VEP_PLUGIN_DIR/Wildtype.pm" "VEP Wildtype plugin"; chk "$VEP_PLUGIN_DIR/Frameshift.pm" "VEP Frameshift plugin"
fi
if wants rna; then
  chk "$STAR_INDEX_DIR/SA" "STAR index"; chk "$SALMON_INDEX_DIR/versionInfo.json" "salmon index"; chk "$TX2GENE" "tx2gene table"
fi
if wants neoantigen; then
  if [[ "$PVAC_IEDB_DIR" == "none" || -z "${PVAC_LOCAL_IEDB_HOST_DIR}" && "$PVAC_IEDB_DIR" != "auto" && "$PVAC_IEDB_DIR" != /* ]]; then
    note "pVACseq will use the IEDB REST API (needs internet, slower)"
  fi
  [[ "${PVAC_ALGORITHMS}" == *NetMHCpan* ]] && note "NetMHCpan runs via a local IEDB install if present in the image, otherwise the IEDB web API. A local DTU-licensed install can be mounted via PVAC_LOCAL_IEDB_HOST_DIR."
fi

info "preflight finished: ${ERRORS} error(s), ${WARNINGS} warning(s)"
(( ERRORS == 0 )) || die "preflight failed — fix the errors above"
