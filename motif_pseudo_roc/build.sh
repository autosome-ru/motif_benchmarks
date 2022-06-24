#!/usr/bin/env bash
set -e -u -o pipefail
BUILD_VERSION=$1
IMAGE_DIR=$(dirname "${BASH_SOURCE[0]}")
docker build \
       -t vorontsovie/motif_pseudo_roc:latest \
       -t vorontsovie/motif_pseudo_roc:"${BUILD_VERSION}" \
       --build-arg BUILD_DATE="$(date -u +'%Y-%m-%dT%H:%M:%SZ')" \
       --build-arg BUILD_VERSION="${BUILD_VERSION}" \
       -f "${IMAGE_DIR}/Dockerfile" \
       "${IMAGE_DIR}"
