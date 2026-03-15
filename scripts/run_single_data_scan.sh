#!/usr/bin/env bash
set -euo pipefail

# Run the uncontrolled ERD data-variation scan from the repository root.
# Extra CLI arguments are passed through to the Python scan script, so you can
# select variants, switch presets, or override the base config as needed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

uv run python erd_fipy/scripts/scan_single_data_variants.py --preset report --tag data_scan "$@"
