# Upstream import — Yul Leshem's analysis code

**Source:** `YLeshem18/nmd_lungcells_2026` (git remote `yul`), branch `main`,
commit `78f1a8c69f8126915dca7ecfeb36f4fb1b353170`, imported 2026-07-20.

**Method:** files only, NOT git history. The two repos have unrelated histories
(`git merge-base` empty), so this is an import + curate, not a merge. To refresh:
re-run `git archive yul/main | tar -x -C code/upstream` after fetching `yul`.

**This is the UPSTREAM half of the pipeline** (Yul owns §1–3 + the DGE/DIE/mashr
infrastructure); the downstream half (Isopair §4, DL model §5, all multipanel
figures, SF25–43) lives in the rest of this repo. See `paper/results_to_code_map.md`
and `REPO_CONSOLIDATION_PLAN.md`.

## Deliberately EXCLUDED from the import
- **`tan_reanalysis/data/*.xlsx`** (4 Tan et al. 2025 supplementary tables) — these
  are another paper's published data, tracked in Yul's repo via **Git LFS**. Not
  imported: (a) redistributing them in a public DOI archive is an open rights question
  (plan W-2), (b) LFS objects import as useless pointer stubs anyway. The file→sheet
  **download manifest is in `tan_reanalysis/README.md`** — download from the Tan paper
  and place in `tan_reanalysis/data/` to run.
- **`.gitignore`, `.gitattributes`** — Yul's repo config (incl. the LFS filter),
  not applicable here.

## Included data
- `data/encode_rbp_roster_vannostrand2020.csv`, `data/gerstberger_2014_rbp_census.csv`
  — RBP rosters for the M8a enrichment analysis (published, small, safe to ship).

## Known follow-ups (from the consolidation plan)
- **D3 path normalization:** these scripts hardcode `/udd/reyle/nmd_lungcells_2026`
  and `.libPaths("/udd/reyle/Rlibs")`, and several `make_*.R` write to a hardcoded
  `HERE`. Must be rewritten to relative paths before the citable snapshot.
- **Data bundle:** the figure scripts read an `nmd_fig_data/` bundle not in either
  repo (deposit workstream D2).
