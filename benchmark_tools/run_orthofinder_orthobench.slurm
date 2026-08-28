#!/usr/bin/env bash
set -euo pipefail

: "${SEARCH:?SEARCH must be diamond or mmseqs}"

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

exec "${SCRIPT_DIR}/run_orthofinder_orthobench.sh" "${SEARCH}"
