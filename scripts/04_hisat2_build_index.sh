#!/usr/bin/env bash
# scripts/04_hisat2_build_index.sh
# Build HISAT2 index for Hordeum vulgare MorexV3 reference
# Usage: bash scripts/04_hisat2_build_index.sh

set -euo pipefail

# Load project config (user-specific paths live here)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
. "${SCRIPT_DIR}/00_config.sh"

mkdir -p "${HISAT2_INDEX_DIR}"

# If you are on an HPC with environment modules, use this.
# Otherwise, comment it out and ensure hisat2-build is on PATH.
module load igmm/apps/HISAT2/2.1.0

# Build index (creates files: MorexV3.1.ht2 ... MorexV3.8.ht2)
hisat2-build "${FASTA}" "${HISAT2_INDEX_BASE}"

echo "HISAT2 index built:"
ls -lh "${HISAT2_INDEX_DIR}" | sed -n '1,200p'
