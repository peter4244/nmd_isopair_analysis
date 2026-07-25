#!/usr/bin/env bash
# =============================================================================
# capture_verification_baseline.sh — regression baseline for the 4-CT-native
# rewrite (REPRODUCIBILITY_PLAN.md §2b.3).
#
# Runs every §1-§3 verification script and stores its normalized stdout. After
# the rewrite, re-run with --check: EVERY number must be identical. Any diff is
# a defect in the rewrite until proven otherwise.
#
#   ./capture_verification_baseline.sh            # write baseline
#   ./capture_verification_baseline.sh --check    # compare against baseline
#
# Normalization strips only volatile lines (absolute paths, timings, package
# load chatter) — never numeric output.
# =============================================================================
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
MODE="${1:-capture}"
BASE="reproducibility/baseline_verification"
OUT="$BASE"; [ "$MODE" = "--check" ] && OUT="tmp/baseline_check"
mkdir -p "$OUT"

SCRIPTS=(
  code/upstream/verify_section3_p1.R
  code/upstream/verify_section3_p1_metric_diagnostics.R
  code/upstream/verify_section3_p2.R
  code/upstream/verify_section3_p3.R
  code/upstream/verify_section3_p4.R
  code/upstream/verify_section3_p5.R
  code/upstream/verify_section3_p2_claim311.R
  code/upstream/verify_section2_claims_210_211.R
  code/upstream/verify_section1_claim113.R
)

# strip volatile content; keep every line that could carry a number
normalize() {
  grep -vE '^\[1\] "/|^Loading required|^ *$|Attaching package|masked from|^purled |^executing |^patched |^applied |^\[runner\]|tidyverse|^✔|^✖|^──|^ℹ|^Warning|^There were|^ *[0-9]+: In |mbcsToSbcs|geom_smooth|select\(\) returned|^saved|^spliced|^committed' \
  | sed -E 's#/Users/[^ ",]*##g; s#[0-9]{4}-[0-9]{2}-[0-9]{2}##g; s#[0-9]+\.[0-9]+ (secs|mins)##g'
}

echo "=== ${MODE} → $OUT ==="
fail=0
for s in "${SCRIPTS[@]}"; do
  n=$(basename "$s" .R)
  printf '%-42s ' "$n"
  if [ ! -f "$s" ]; then echo "MISSING"; fail=1; continue; fi
  Rscript "$s" 2>&1 | normalize > "$OUT/$n.txt"
  lines=$(wc -l < "$OUT/$n.txt" | tr -d ' ')
  if [ "$MODE" = "--check" ]; then
    if diff -q "$BASE/$n.txt" "$OUT/$n.txt" >/dev/null 2>&1; then
      echo "IDENTICAL ($lines lines)"
    else
      echo "*** DIFFERS ***"; fail=1
      diff "$BASE/$n.txt" "$OUT/$n.txt" | head -25 | sed 's/^/      /'
    fi
  else
    echo "captured ($lines lines)"
  fi
done

echo
if [ "$MODE" = "--check" ]; then
  [ $fail -eq 0 ] && echo "✅ ALL IDENTICAL — rewrite is numerically faithful" \
                  || echo "❌ DIFFERENCES FOUND — investigate before proceeding"
else
  echo "Baseline written to $BASE/. Commit it BEFORE editing any analysis code."
fi
exit $fail
