#!/usr/bin/env bash
# Build local bioconda-based images and pull the official ones.
set -Eeuo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
init_step "00_build_images"

command -v docker >/dev/null 2>&1 || die "docker not found in PATH"

for img in "${!LOCAL_IMAGES[@]}"; do
  pkgs="${LOCAL_IMAGES[$img]}"
  if image_present "$img" && (( ! FORCE )); then info "image present: $img"; continue; fi
  info "building $img  (packages: $pkgs)"
  if (( DRY_RUN )); then continue; fi
  CUR_LOG="$LOG_DIR/build_${img//[\/:]/_}.log"
  if ! docker build -f "$DOCKER_CTX_DIR/Dockerfile.conda" \
        --build-arg "BASE=$MICROMAMBA_BASE" --build-arg "PACKAGES=$pkgs" \
        -t "$img" "$DOCKER_CTX_DIR" >"$CUR_LOG" 2>&1; then
    die "docker build failed for $img (see $CUR_LOG). A package pin may have moved: edit LOCAL_IMAGES in config/config.sh"
  fi
  CUR_LOG=""
done

for img in "$IMG_GATK" "$IMG_VEP" "$IMG_PVAC" "$IMG_BAMRC"; do
  if (( DRY_RUN )); then info "[dry-run] would pull $img"; continue; fi
  if image_present "$img" && (( ! FORCE )); then info "image present: $img"; continue; fi
  info "pulling $img (pVACtools is ~10 GB compressed; this takes a while)"
  docker pull "$img" || die "Cannot pull $img. Edit its tag in config/config.sh (IMG_*)."
done
info "images ready"
