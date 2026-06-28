#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba"
VARIANTS_FILE="$PKG_DIR/a1_gpiba_all_recommended_variants.txt"   bash "$PKG_DIR/run_a1_gpiba_p0_md.sh"
