#!/usr/bin/env bash
set -euo pipefail

# Remove leftover temp artifacts from tests

# 1) mkdtemp-created baseline dirs under /tmp
find /tmp -maxdepth 1 -type d -name 'solo_baseline_test.*' -print -exec rm -rf {} + || true

# 2) Repo-local baseline capture outputs (if any)
REPO_ROOT=$(cd "$(dirname "$0")/../.." && pwd)
rm -rf "$REPO_ROOT/new/testing/outputs/baseline_skip_mode" || true

# 3) Any stray _STARtmp under repo testing outputs
find "$REPO_ROOT/new/testing/outputs" -type d -name '_STARtmp' -print -exec rm -rf {} + 2>/dev/null || true

echo "Cleanup complete."

