#!/usr/bin/env bash
# =============================================================================
#  run.sh — entry point of the melanoma neoantigen / mRNA design pipeline
#  RESEARCH USE ONLY.   See README.md.   Run ./run.sh --help
# =============================================================================
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$ROOT/scripts"
# shellcheck disable=SC1091
source "$ROOT/config/config.sh"

usage() {
cat <<'EOF'
Usage: ./run.sh [options] [STAGE ...]

Stages (default: all)
  setup        build/pull Docker images + download & index references (once, ~2-4 h, ~100 GB)
  qc           01 FastQC + 02 fastp
  align        03 BWA alignment + 04 MarkDuplicates/coverage
  variants     05 Mutect2 + 06 VEP annotation
  rna          07 STAR + 08 Salmon
  hla          09 OptiType (class I, from normal exome)
  neoantigen   10 pVACseq (NetMHCpan...) + 11 ranking + 12 mRNA design
  report       13 MultiQC + HTML report
  all          preflight + qc align variants rna hla neoantigen report (in numeric step order)
  status       show which units are complete
  cleanup      free disk space (see --level)

Options
  --step N[,N..]       run only numbered step(s), e.g. --step 05,06
  --from N --to N      run a numeric range of steps, e.g. --from 03 --to 06
  --patient ID         restrict to one patient from the samplesheet
  --samplesheet FILE   alternative samplesheet CSV
  --force              ignore resume markers and recompute selected stages/steps
  --dry-run            print the Docker commands without running anything
  --keep-intermediate  do not delete trimmed FASTQ / merged BAMs / shards
  --skip-preflight     do not run 00_preflight first
  --components LIST    for 'setup': images genome gatk gencode bwa star salmon vep (default all)
  --level LEVEL        for 'cleanup': tmp | intermediate | docker | all  (default tmp)
  -y, --yes            do not ask for confirmation in cleanup
  -h, --help

Examples
  ./run.sh setup                      # once
  ./run.sh                            # whole pipeline, resumes automatically
  ./run.sh variants rna               # two stages
  ./run.sh --step 10,11 --patient Patient01 --force
  ./run.sh cleanup --level intermediate
EOF
}

declare -A STAGE_STEPS=(
  [qc]="01 02" [align]="03 04" [variants]="05 06" [rna]="07 08" [hla]="09"
  [neoantigen]="10 11 12" [report]="13"
)
STAGE_ORDER=(qc align variants rna hla neoantigen report)
declare -A STEP_FILE=(
  [01]=01_fastqc [02]=02_fastp [03]=03_bwa_mem2 [04]=04_markdup [05]=05_mutect2 [06]=06_annotation
  [07]=07_star [08]=08_salmon [09]=09_hla [10]=10_netmhcpan [11]=11_ranking [12]=12_mrna_design [13]=13_report
)

stages=(); steps=(); from=""; to=""; skip_pre=0; components="all"; level="tmp"; yes=0
export FORCE=0 DRY_RUN=0 PATIENT_FILTER="" KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
while (( $# )); do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --force) FORCE=1 ;;
    --dry-run) DRY_RUN=1 ;;
    --keep-intermediate) KEEP_INTERMEDIATE=1 ;;
    --skip-preflight) skip_pre=1 ;;
    --patient) PATIENT_FILTER="$2"; shift ;;
    --samplesheet) export SAMPLESHEET="$2"; shift ;;
    --step) IFS=',' read -ra s <<< "$2"; for x in "${s[@]}"; do steps+=("$(printf '%02d' "$((10#$x))")"); done; shift ;;
    --from) from="$(printf '%02d' "$((10#$2))")"; shift ;;
    --to) to="$(printf '%02d' "$((10#$2))")"; shift ;;
    --components) components="$2"; shift ;;
    --level) level="$2"; shift ;;
    -y|--yes) yes=1 ;;
    -*) echo "unknown option: $1" >&2; usage; exit 2 ;;
    *) stages+=("$1") ;;
  esac
  shift
done
export FORCE DRY_RUN PATIENT_FILTER KEEP_INTERMEDIATE

mkdir -p "$LOG_DIR" "$RESULTS_DIR" "$TMP_DIR"
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"; export MELPIPE_RUN_ID="$RUN_ID"
MAIN_LOG="$LOG_DIR/run_${RUN_ID}.log"
exec > >(tee -a "$MAIN_LOG") 2>&1
export MELPIPE_LOGGING=1
say() { printf '[%s] [run.sh] %s\n' "$(date '+%F %T')" "$*"; }

# ---- single-instance lock ----------------------------------------------------
exec 9>"$RESULTS_DIR/.run.lock"
flock -n 9 || { say "another pipeline run holds the lock ($RESULTS_DIR/.run.lock)"; exit 1; }

# ---- stop containers when interrupted -----------------------------------------
on_signal() {
  say "interrupted — stopping containers started by this run"
  docker ps -q --filter "label=melpipe_run=${RUN_ID}" 2>/dev/null | xargs -r docker kill >/dev/null 2>&1 || true
  exit 130
}
trap on_signal INT TERM

# ---- special commands ----------------------------------------------------------
[[ ${#stages[@]} -eq 0 && ${#steps[@]} -eq 0 && -z "$from$to" ]] && stages=(all)

run_script() {   # run_script FILE   (returns the script's exit code)
  local f="$1" t0 rc
  t0=$(date +%s)
  say "▶ $f"
  set +e
  bash "$SCRIPTS/$f.sh" 2>&1 | tee -a "$LOG_DIR/${f}.log"
  rc=${PIPESTATUS[0]}
  set -e
  SUMMARY+=("$(printf '%-16s %-8s %s' "$f" "$( ((rc==0)) && echo OK || echo FAILED )" "$(fmt_dur $(( $(date +%s) - t0 )))")")
  return "$rc"
}
fmt_dur() { printf '%dh%02dm%02ds' $(($1/3600)) $(($1%3600/60)) $(($1%60)); }
SUMMARY=()
print_summary() {
  say "================ summary ================"
  printf '%s\n' "${SUMMARY[@]:-}" | sed 's/^/    /'
  say "main log: $MAIN_LOG"
}

for s in "${stages[@]:-}"; do
  [[ -n "$s" ]] || continue
  case "$s" in
    setup)
      SETUP_COMPONENTS="${components//,/ }"; export SETUP_COMPONENTS
      if [[ "$components" == "all" || " ${SETUP_COMPONENTS} " == *" images "* ]]; then run_script 00_build_images || { print_summary; exit 1; }; fi
      if [[ "$components" != "images" ]]; then run_script 00_setup_reference || { print_summary; exit 1; }; fi
      print_summary; exit 0 ;;
    cleanup)
      if (( yes )); then export CLEAN_YES=1; fi
      bash "$SCRIPTS/cleanup.sh" "$level"; exit 0 ;;
    status)
      say "markers in $MARK_DIR:"
      for m in "$MARK_DIR"/*.done; do
        if [[ -e "$m" ]]; then printf '  %-48s %s\n' "$(basename "$m" .done)" "$(cat "$m")"; fi
      done
      say "disk: results=$(du -sh "$RESULTS_DIR" 2>/dev/null | cut -f1) reference=$(du -sh "$REF_DIR" 2>/dev/null | cut -f1) free=$(df -Pk "$PROJECT_ROOT" | awk 'NR==2{printf "%d GB", $4/1024/1024}')"
      exit 0 ;;
    all|qc|align|variants|rna|hla|neoantigen|report) ;;
    *) say "unknown stage '$s'"; usage; exit 2 ;;
  esac
done

# ---- resolve the ordered list of numbered steps ------------------------------------
want=()
if (( ${#steps[@]} )); then want=("${steps[@]}")
elif [[ -n "$from$to" ]]; then
  for n in 01 02 03 04 05 06 07 08 09 10 11 12 13; do
    [[ -n "$from" && "$n" < "$from" ]] && continue
    [[ -n "$to" && "$n" > "$to" ]] && continue
    want+=("$n")
  done
else
  chosen=()
  for s in "${stages[@]}"; do
    if [[ "$s" == all ]]; then chosen=("${STAGE_ORDER[@]}"); break; fi
    chosen+=("$s")
  done
  for n in 01 02 03 04 05 06 07 08 09 10 11 12 13; do
    for s in "${chosen[@]}"; do
      [[ " ${STAGE_STEPS[$s]} " == *" $n "* ]] && { want+=("$n"); break; }
    done
  done
fi
(( ${#want[@]} )) || { say "nothing selected"; exit 2; }
for n in "${want[@]}"; do [[ -n "${STEP_FILE[$n]:-}" ]] || { say "invalid step '$n'"; exit 2; }; done

# stage names needed for preflight
pf=()
for n in "${want[@]}"; do
  for s in "${STAGE_ORDER[@]}"; do [[ " ${STAGE_STEPS[$s]} " == *" $n "* ]] && pf+=("$s"); done
done
export PREFLIGHT_STAGES="$(printf '%s\n' "${pf[@]}" | awk '!s[$0]++' | tr '\n' ' ')"

say "run id $RUN_ID | steps: ${want[*]} | force=$FORCE dry_run=$DRY_RUN patient='${PATIENT_FILTER:-all}'"
if (( ! skip_pre )); then
  run_script 00_preflight || { print_summary; say "preflight failed"; exit 1; }
fi
for n in "${want[@]}"; do
  if ! run_script "${STEP_FILE[$n]}"; then
    print_summary
    say "step $n FAILED. Fix the cause and re-run the same command: finished units are skipped automatically."
    exit 1
  fi
done
print_summary
say "all requested steps finished"
