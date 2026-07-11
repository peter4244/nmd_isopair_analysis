# Supplemental-Figure legend QA checklist

**Run this on EVERY SF legend before rebuilding the docx.** These are the corrections Pete has had to flag repeatedly — every item below has bitten a real legend. Mechanical items (§E) are enforced by `figures/lib/lint_sf_legends.py`; run it first, then eyeball the judgment items.

> Fast path: `python3 figures/lib/lint_sf_legends.py` → fix every hit → then walk §A–§D for the figure(s) you touched.

---

## A. Numbers in the legend MUST match the figure — the #1 recurring bug
The legend and the rendered figure are authored/generated separately, so they drift. **Every count, median, percentage, and p-value stated in prose must equal what the figure/data show.**

- [ ] **Sample sizes (`n = …`)** in the legend == the n's the figure renders (facet headers, x-tick `(n=…)`, panel annotations). This is the one that keeps breaking: SF39 legend said 1,016/95/1,107 while the figure showed 794/60/851; same class caught in SF40 (255/30/276 → 202/20/219) and SF32 (72/118/190 → 54/82/136). **Pull the number from the figure's own data TSV, never from memory or an older legend.**
- [ ] Model-subset / scope-specific figures (SF39/40, GENCODE-restricted SF32/35, occult SF34) cite *their own* subset counts — not the pop_BC headline numbers. A grep for the headline n's won't catch a stale subset n; check each explicitly.
- [ ] Medians / IQRs / percentages / OR / p-values in prose match the panel annotations and the descriptives TSV.
- [ ] Denominators are right: "% of X" states the correct X, and the population compared is the intended one.

## B. Terminology & symbols (project conventions — enforced by lint)
- [ ] **"exon junction complex" → "EJC"** everywhere in prose, AND `EJC, exon junction complex` appears in that legend's **Abbreviations**. (Recurring.)
- [ ] **Start/stop codons: "AUG" not "ATG"** in RNA-side display text, axis labels, and legends (data-side column names may stay "ATG"). SF33/34 shipped with "ATG" on an axis. "start codon"/"stop codon" are **never hyphenated**, even as modifiers; prefer "start codon" over "AUG" in user-facing *titles*.
- [ ] **Significance stars must be defined** whenever used: e.g. `*p < 0.05, **p < 10⁻³, ***p < 10⁻⁴; n.s. otherwise`. A legend that says "significance stars" or shows `*/**/***` without the key is incomplete. (SF29 recurred.) **The key must define every level the figure actually renders** — if a panel draws `****` (p<10⁻¹⁰, e.g. SF32/35 with p≈10⁻⁴⁹), a `*/**/***`-only key is incomplete; use the canonical 4-level `p_to_stars` ladder.
- [ ] **SF-legend star keys must ESCAPE the asterisks** (`\*p < 0.05, \*\*p < 10⁻³, \*\*\*p < 10⁻⁴, \*\*\*\*p < 10⁻¹⁰`). SF legends compile into the **markdown** SF-docx (`build_supplemental_figures_docx.js`), which otherwise eats bare `**`/`***`/`****` as bold and drops the leading stars (SF29 shipped mangled 2026-07-11). Now enforced by `lint_sf_legends.py`. **Exception:** the multipanel Fig 3/4/5 composite legends stay **unescaped** — they feed the manuscript Google Doc as plain text, where a literal `*` is what you want.
- [ ] **"NMD susceptible"** — no hyphen — in prose (R column `nmd_responsive` unchanged).
- [ ] **US spelling**: color / center / gray / analyze (no coloured/centred/grey/analyse).
- [ ] **No project jargon** in paper-facing legend text: no `pop_BC`, `pair-set`, `profiles_c2`, chunk/Rmd names, category codes. Use accepted terms: "main open reading frame" (not priority ORF), "occult PTC" (not hidden PTC), "TransDecoder2".

## C. Style & scope (see `~/.claude/memory/figures_style_publications.md`; read it BEFORE drafting)
- [ ] **Title Case** for the `SF## | …` title, preserving domain acronyms (CDS, NMD, PTC, |SHAP|, ORF; `mashr` stays lowercase).
- [ ] **No back-references to manuscript sections** ("as described in §5") — legends are standalone.
- [ ] **No narration of the visual shape** ("the violin is bimodal"), and **no mechanism/biology interpretation**. State what is plotted, not what it means.
- [ ] Panel-by-panel (A/B/C…) coverage; every panel referenced; line/marker keys defined (dashed = median, etc.) and **matching what the figure actually draws** — SF26 legend named a "dotted 50% line" the figure no longer had.
- [ ] Use Pete's verbatim prose when he supplies canonical legend text — sync exactly, don't paraphrase.

## D. Basis / anchor documentation — REQUIRED, but only the methodological choice
Legends must document *which coordinate/subset/basis* is shown (this is the one kind of "methods" a legend should carry):
- [ ] **Stop/PTC/CDS anchor source** stated when relevant: ref-AUG-traced (unbiased) vs TD2-called (skips PTCs). Any post-stop / PTC / CDS measure must name its anchor; prefer ref-AUG-anchored. (See `feedback_nmd_paper_check_cds_source`.)
- [ ] **Expression basis is 4-CT.** After the 2026-07 re-scope, the whole analysis is AT/DD/FB/MV only. Legends say **"the four cell types"**, never "all sequenced libraries", "4-cell-type across AT2/LAE/FB/MV" as a basis claim that's actually 6-CT, or anything implying DD_ALI/DO. No 6-CT trace anywhere.
- [ ] Denominator basis stated correctly (e.g. reference share is `ref / Σ all isoforms`, i.e. **total gene expression**, not "non-NMD expression").

## E. Mechanical pre-flight — `figures/lib/lint_sf_legends.py`
Flags, across all `*legend*.md`: `exon junction complex` without an EJC abbreviation entry; `ATG` in display text; `significance star`/`*`-usage with no threshold key; `all sequenced libraries` / `DD_ALI` / `DO` basis leaks; `non-NMD expression` as a denominator; sentence-case SF titles; `§`/section back-refs; British spellings. **A clean lint is necessary, not sufficient — §A number-matching is still manual (it needs the figure data).**

---
*Provenance: consolidated 2026-07-11 from repeat corrections (SF39/40/32 n-mismatch, SF26 denominator/basis/line, SF33/34 ATG, SF29 stars, EJC). Related memories: `feedback_nmd_sf_legend_style`, `feedback_figure_legend_read_style_first`, `feedback_nmd_paper_check_cds_source`, `feedback_nmd_sf_title_case`, `feedback_use_user_provided_prose_verbatim`.*
