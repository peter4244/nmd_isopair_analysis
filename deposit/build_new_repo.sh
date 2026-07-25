#!/usr/bin/env bash
# =============================================================================
# build_new_repo.sh — assemble nmd_lung_longread_2026 from the manifest.
#
# Copies only git-TRACKED files from the two source repos, so nothing untracked
# or ignored leaks in. Idempotent: wipes and rebuilds the target each run.
# Layout and exclusions follow NEW_REPO_MANIFEST.md.
# =============================================================================
set -euo pipefail
SRC="$HOME/claude_projects/nmd"
MODEL="$HOME/claude_projects/NMD_orf_model_v5_4ct"
DEST="$HOME/claude_projects/nmd_lung_longread_2026"

rm -rf "$DEST"; mkdir -p "$DEST"
mkdir -p "$DEST"/{analysis/{upstream,isopair,predictor_comparison},model,figures,verification/baseline,metadata/{pheno,reference},docs,config}

# copy a list of tracked paths (relative to a repo) preserving structure below a prefix
copy_tracked() { # $1=repo $2=pathspec $3=strip-prefix $4=dest
  local repo="$1" spec="$2" strip="$3" dest="$4" n=0
  cd "$repo"
  while IFS= read -r f; do
    [ -f "$f" ] || continue
    local rel="${f#$strip}"
    mkdir -p "$dest/$(dirname "$rel")"
    cp "$f" "$dest/$rel"; n=$((n+1))
  done < <(git ls-files "$spec")
  echo "$n"
}

echo "=== analysis/upstream (§1–§3) ==="
n=$(copy_tracked "$SRC" "code/upstream" "code/upstream/" "$DEST/analysis/upstream"); echo "  $n files"
# verification scripts and reference rosters relocate out of upstream/
mkdir -p "$DEST/verification"
mv "$DEST/analysis/upstream"/verify_*.R "$DEST/verification/" 2>/dev/null || true
if [ -d "$DEST/analysis/upstream/data" ]; then rm -rf "$DEST/analysis/upstream/data"; fi
echo "  -> moved $(ls "$DEST/verification" | wc -l | tr -d ' ') verify_*.R to verification/"

echo "=== analysis/isopair (§4) ==="
n=$(copy_tracked "$SRC" "results/isoform_transitions/Version_6.0" \
      "results/isoform_transitions/Version_6.0/" "$DEST/analysis/isopair"); echo "  $n files"
rm -rf "$DEST/analysis/isopair/results"          # rmats/PTC data products — excluded
echo "  -> removed results/ptc (excluded per manifest)"

echo "=== analysis/predictor_comparison (§5 benchmark) ==="
n=$(copy_tracked "$SRC" "code/nmd_predictor_comparison" "code/nmd_predictor_comparison/" \
      "$DEST/analysis/predictor_comparison"); echo "  $n files"

echo "=== model (§5) ==="
cd "$MODEL"
n=0; while IFS= read -r f; do
  case "$f" in superseded/*) continue;; esac
  [ -f "$f" ] || continue
  mkdir -p "$DEST/model/$(dirname "$f")"; cp "$f" "$DEST/model/$f"; n=$((n+1))
done < <(git ls-files)
echo "  $n files (superseded/ excluded)"

echo "=== figures ==="
n=$(copy_tracked "$SRC" "figures" "figures/" "$DEST/figures"); echo "  $n files"
[ -d "$DEST/figures/SupplementalFigures" ] && mv "$DEST/figures/SupplementalFigures" "$DEST/figures/supplemental"

echo "=== verification (verifiers + baseline) ==="
n=$(copy_tracked "$SRC" "reproducibility" "reproducibility/" "$DEST/verification"); echo "  $n files"

echo "=== metadata ==="
n=$(copy_tracked "$SRC" "pheno" "pheno/" "$DEST/metadata/pheno"); echo "  pheno: $n"
cp "$SRC"/code/upstream/data/*.csv "$DEST/metadata/reference/" 2>/dev/null || true
echo "  reference rosters: $(ls "$DEST/metadata/reference" 2>/dev/null | wc -l | tr -d ' ')"

echo "=== docs ==="
cp "$SRC/paper/results_to_code_map.md" "$DEST/docs/" 2>/dev/null || true
cp "$SRC/ENVIRONMENT.md" "$DEST/" 2>/dev/null || true

echo
echo "=== TOTAL: $(find "$DEST" -type f | wc -l | tr -d ' ') files ==="
find "$DEST" -type f | sed "s#$DEST/##" | awk -F/ '{print $1"/"($2 ~ /\./ ? "" : $2)}' | sort | uniq -c | sort -rn
