#!/usr/bin/env bash
# Storage clean-up.   scripts/cleanup.sh [tmp|intermediate|docker|all] [--yes]
#   tmp           empty data/tmp (sort/STAR scratch) and stale *.partial files
#   intermediate  + trimmed FASTQ, merged/lane BAMs, Mutect2 work dirs, pVACseq prep dirs
#   docker        + dangling Docker layers/build cache (images are kept)
#   all           everything above
# Final deliverables (md.bam, VCFs, RNA BAM, tables, reports) are never touched.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
STEP_NAME="cleanup"; mkdir -p "$LOG_DIR"
level="${1:-tmp}"; yes=0; [[ "${2:-}" == "--yes" || "${CLEAN_YES:-0}" == 1 ]] && yes=1

report() { info "disk: results=$(du -sh "$RESULTS_DIR" 2>/dev/null | cut -f1) reference=$(du -sh "$REF_DIR" 2>/dev/null | cut -f1) tmp=$(du -sh "$TMP_DIR" 2>/dev/null | cut -f1) free=$(free_gb "$PROJECT_ROOT") GB"; }
report
targets=()
case "$level" in
  tmp|intermediate|docker|all) ;;
  *) die "usage: cleanup.sh [tmp|intermediate|docker|all] [--yes]" ;;
esac
targets+=("$TMP_DIR"/*)
while IFS= read -r f; do targets+=("$f"); done < <(find "$RESULTS_DIR" -name '*.partial*' 2>/dev/null)
if [[ "$level" == intermediate || "$level" == all ]]; then
  for d in "$RESULTS_DIR"/trimmed/* ; do targets+=("$d"); done
  while IFS= read -r f; do targets+=("$f"); done < <(find "$RESULTS_DIR/alignment" \( -name '*.merged.bam' -o -name '*.sorted.bam' \) 2>/dev/null)
  while IFS= read -r f; do targets+=("$f"); done < <(find "$RESULTS_DIR/variants" -type d -name mutect2_work 2>/dev/null)
  while IFS= read -r f; do targets+=("$f"); done < <(find "$RESULTS_DIR/neoantigen" -type d -name prep 2>/dev/null)
fi
existing=(); for t in "${targets[@]}"; do [[ -e "$t" ]] && existing+=("$t"); done
if (( ${#existing[@]} )); then
  printf '  would remove: %s\n' "${existing[@]}" >&2
  if (( ! yes )); then read -r -p "Delete the items above? [y/N] " a; [[ "$a" == y || "$a" == Y ]] || { info "aborted"; exit 0; }; fi
  rm -rf -- "${existing[@]}"
else info "nothing to remove"; fi
if [[ "$level" == docker || "$level" == all ]]; then docker system prune -f >/dev/null 2>&1 && info "docker prune done" || warn "docker prune failed"; fi
mkdir -p "$TMP_DIR"; report
