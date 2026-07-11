# Session state (before context compaction)

## Where we are

Working through Pete's collected defect list (13 items) while validator architecture planning happens in parallel.

## Architectural planning (in-flight)

- **Plan v1** at `/private/tmp/claude-502/-Users-petecastaldi/66ba122b-82af-4cd2-b3cf-02f55870509d/scratchpad/validator_plan_v1.md`. Recommends Architecture A (typed element enumerator + pairwise rule registry + RenderContext value object), incremental 6-phase migration over ~7 days. Pete approved the direction but asked for a critique.
- **Critic agent (agentId a856808906343d7f4)** running in background — adversarial review of plan v1. Notification will fire when done. Pete's plan: review v1 + critique together, then implement.
- **Fork Pete still needs to decide** before Phase 3: (a) sidecar `WeakKeyDictionary` for RenderContext-on-fig (plan-recommended), (b) attribute `fig._nmd_render_context`.

## Defect list (13 items)

Numbered per collection order. **Status codes:** ✅ done, 🔄 in progress, ⏸ blocked on Pete input, ⏳ pending.

1. **✅ SF27 legend** — removed "dotted bars, 25th and 75th percentiles" (this session)
2. **✅ SF29** — 50% line + text removed (2a); "GAIN" → "Gain" (2b); legend moved to upper center. Committed in `7493387`.
3. **✅ Build script** referenced legacy dirs — fixed in commit `2f7d031`; docx now has correct separate SF32/33/34/35
4. **⏳ Delete legacy dirs** — needs `reproducibility/verify_cross_check_new_rmd_vs_figures.R` paths updated first (lines 26-28, 147, 154, 161 reference SUPP_CDS = SF32-SF35 legacy dir, SUPP_TD2 = SF33-SF34 legacy dir)
5. **🔄 SF32 legend** — `\*p < 0.05, \*\*p < 10⁻³, \*\*\*p < 10⁻⁴` — asterisks are backslash-escaped in markdown source. Also SF35 has the same pattern. Fix: change `\*` → `*` in both files. Consider docx-build-script implications for other markdown special chars.
6. **⏳ SF33 Panel C** — xtick `0`/`−500` overlap at origin. Fix: set explicit sparser xticks (like SF30 had), likely `[-500, 0, 500, 1000, 1500, 2000, 2500]`.
7. **⏳ SF34 Panel C** — same overlap. Same fix (companion figure — always sync).
8. **⏳ SF37 legend** — "Abbreviations" on separate line. Fix: put inline like other SFs.
9. **⏳ Repo-wide** — `N %` (space between number and %) → `N%`. Grep `figures/SupplementalFigures/**/*_legend.md`.
10. **⏳ SF38 legend prose** — too long AND contains interpretive commentary. Violates the `feedback_nmd_sf_legend_style.md` rule (legends must not interpret biology/mechanism, just describe what's shown). Effectively subsumed by #12 — extracting the brief docx version replaces it with the non-interpretive short form. Do #12 first; if any interpretive prose survives after extraction, strip it manually.
11. **⏸ SF38 in-plot legend overlaps content** — same class as #2a (validator gap). Handled by B3 validator (assert_legend_clear) but SF38 still fails with ~5000 px overlap; part of the 4 residual B3 failures.
12. **⏳ SF38-SF42 legends: verbose style** — extract legends from **base `nmd_supplemental_figures_sf24_sf42.docx`** (Pete confirmed 2026-07-10) and replace the on-disk `_legend.md` files for SF38, SF39, SF40, SF41, SF42.
13. **⏳ SF29 xtick crowding** — "Missing Internal" (2-line) and "IR diff 3′" touch. Validator gap: `assert_tick_labels_disjoint` needs `min_gap_px=~8` parameter, not just no-overlap. Would also close #6/#7 gaps.

## Validator work still queued (from B roadmap)

- **B1** (legend-vs-other-text): NOT started. Deferred by architectural refactor — Architecture A's element-graph should handle this cleanly. Consider skipping B1 as standalone and doing it as first new check under new architecture.
- **B2** (tick_labels_disjoint min_gap_px param): NOT started. Small fix that could be done now OR deferred to post-refactor.
- **B3** (legend_clear): ✅ landed in commit `7493387`. 4 residual figure failures (SF36, SF38, SF40, SF41) — small overlaps 1900-5000 px, density-tail-in-corner. Pete's guidance: try external legend or per-figure `overlap_tolerance_px` override.

## B3 residuals — 4 SFs still fail render

SF36, SF38, SF40, SF41. Each has 1900-5000 px legend-vs-data overlap where a density curve tail reaches every corner. On the roadmap for after list items.

## Session-commit history (relevant to this work)

Recent commits in `~/claude_projects/nmd`:
- `7493387` Validator B3: assert_legend_clear + fix SF29/SF31/SF42
- `2f7d031` docx: point build script at separate SF32/33/34/35 dirs
- `4069984` Rebuild docx with 12 fixed SFs
- `0380afc` Backfill render_and_validate into 11 remaining SFs
- `d06560c` SSOT-owned render_and_validate + 6 SFs
- `5a0ffab` SF42: BODY_FS + assert_docx_readable
- `7550824` SF42: 2-row layout
- `dd121d9` SF40: re-render
- `e7cb4ce` SF33/SF34: Density y-label
- `2bf6ab3` SF32/SF35 companion

## Non-obvious workflow notes

- **Working tree has untracked files** that are Pete's (paper/*.docx PJCcopy, node_modules, isocall_pipeline.channing-backup) — do not touch these.
- **Every git push goes to public GitHub** + Channing GitLab dual-push per CLAUDE.md — review before committing.
- **`~/miniforge3/bin/python`** is the Python interpreter with matplotlib/pandas/seaborn — the base python3 lacks pandas.
- **Cross-project utils** at `~/.claude/utils/`: `figure_template.py` (starter), `figure_lint.py` (checks no-fontsize-literals + render_and_validate present).
- **Workflow doc** at `~/.claude/memory/figures_workflow_publications.md` — prescribes `render_and_validate` as canonical Phase 4.

## Immediate next steps if continuing

1. Complete legend text fixes: #5 (SF32/SF35 backslash), #8 (SF37 Abbrev), #9 (N% grep). Each is a 1-file edit.
2. Complete tick fixes: #6/#7 (SF33/SF34 Panel C xticks). Requires render + gate check.
3. #4 (delete legacy dirs) — needs reproducibility script path update first.
4. When critic agent notification arrives: review both plans together, decide on architecture direction, decide fork (a) vs (b).
5. Ask Pete about #10 and #12 clarifications when convenient.
