# isocall_dge/mashr/

Mashr-based differential isoform expression (DIE) results from the isocall pipeline.

## Files (all gitignored)

- `nmd_mashr_die_{at,dd,ddali,doali,fb,mv}_2026.3.10.csv` — per-cell-type mashr posterior estimates
- `mashr_isoform_model_2026.3.10.rds` — fitted mashr model (~75 MB)

## Cell-type scope

Mashr pipeline includes `doali` (note: aligned with the project-wide `DO → DO_ALI` rename); the `do` label used by the limma side is not present here.

## Why gitignored

CSVs are ~18–19 MB each; the model `.rds` is ~75 MB. Regenerable from the limma output via the mashr fit pipeline, so kept out of git history (same pattern as `longread_dge/*.csv` and `shortread_dge/*.csv`). See `~/.claude/projects/-Users-petecastaldi-claude-projects-nmd/memory/mashr_pipeline.md` for pipeline details.
