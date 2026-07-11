#!/usr/bin/env python3
"""Lint Supplemental-Figure legend .md files for the recurring correction classes.

Enforces the mechanical items of figures/SupplementalFigures/LEGEND_QA_CHECKLIST.md.
Pure stdlib (runs under any python3). It CANNOT check §A number-matching — that
needs the figure data — so a clean run is necessary, not sufficient.

Usage:
    python3 figures/lib/lint_sf_legends.py [path ...]
Default path: figures/SupplementalFigures. Exits 1 if any ERROR-level flag fires.
"""
import re
import sys
from pathlib import Path

SMALL = {"of", "the", "and", "in", "by", "vs", "per", "to", "a", "an", "for",
         "across", "from", "on", "at", "with", "the", "as"}
ACRONYMS_OK = True  # acronyms (CDS, NMD, PTC, ORF, EJC, |SHAP|) are already caps

def check(path: Path):
    text = path.read_text(encoding="utf-8")
    low = text.lower()
    errs, warns = [], []

    # B: exon junction complex without an EJC abbreviation entry
    if "exon junction complex" in low:
        if not re.search(r"EJC\s*,\s*exon junction complex", text, re.I):
            errs.append('"exon junction complex" present but no "EJC, exon junction '
                        'complex" in Abbreviations (substitute EJC + add to Abbreviations)')

    # B: ATG in display text (data-side column names are rare in legends; flag all)
    for m in re.finditer(r"\bATG\b", text):
        errs.append(f'"ATG" in legend text (use "AUG" for RNA-side display) '
                    f'@ …{text[max(0,m.start()-25):m.start()+8]}…')

    # B: significance stars used but no threshold key
    uses_stars = bool(re.search(r"significance star", low)) or \
                 bool(re.search(r"(?<![\w])\\?\*{1,3}\s*p\s*[<=]", text)) or \
                 " n.s." in low
    has_key = bool(re.search(r"\\?\*+\s*p\s*[<=]\s*10|\\?\*+\s*p\s*[<=]\s*0?\.\d", text))
    if uses_stars and not has_key:
        errs.append("significance stars / n.s. referenced but no threshold key "
                    "(add e.g. *p<0.05, **p<10⁻³, ***p<10⁻⁴; n.s. otherwise)")
    # B: star-key thresholds present but the multi-asterisks are UNESCAPED -> the
    #    SF-docx markdown build (build_supplemental_figures_docx.js) eats them as
    #    bold/italic markers (SF29 recurred 2026-07-11: '****p<10⁻¹⁰, ***p<10⁻⁴'
    #    rendered '**p<10⁻⁴**' bold with leading stars dropped). A properly escaped
    #    key ('\*\*\*\*p') has a backslash between every asterisk, so no run of >=2
    #    bare consecutive '*' survives -- flag any such run before 'p<'/'p='.
    #    Scope: SF legends only -- they compile into the markdown SF-docx. The
    #    multipanel Fig3/4/5 composite legends feed the manuscript Google Doc as
    #    plain text, where literal '*' is correct, so they are exempt.
    if has_key and "SupplementalFigures" in str(path) \
            and re.search(r"\*{2,}\s*p\s*[<=]", text):
        errs.append("star-key thresholds use UNESCAPED multi-asterisks "
                    "(e.g. '****p < 10⁻¹⁰') -- the SF-docx markdown build renders "
                    "them bold; escape as '\\*\\*\\*\\*p' like SF32/SF35.")

    # D: 6-CT basis leaks (post-2026-07 everything is 4-CT).
    # Cell-type CODES are uppercase → case-SENSITIVE (else \bDO\b matches the verb "do").
    for pat, msg in [
        (r"\bDD_ALI\b", "DD_ALI reference (out-of-scope cell type; no 6-CT trace)"),
        (r"\bDO_ALI\b", "DO_ALI reference (out-of-scope cell type)"),
        (r"\bDO\b(?![I_])", "'DO' cell-type code (out-of-scope) — verify not a stray match"),
    ]:
        for m in re.finditer(pat, text):  # no re.I
            errs.append(f'{msg} @ …{text[max(0,m.start()-20):m.start()+15]}…')
    # Basis WORDING is case-insensitive.
    for pat, msg in [
        (r"all sequenced libraries", '"all sequenced libraries" basis wording (say "the four cell types")'),
        (r"\b(six|6)[ -]cell", "six-cell-type basis wording (analysis is 4-CT)"),
    ]:
        for m in re.finditer(pat, text, re.I):
            errs.append(f'{msg} @ …{text[max(0,m.start()-20):m.start()+15]}…')

    # D: non-NMD denominator mislabel
    if re.search(r"non-?NMD expression", text, re.I):
        errs.append('"non-NMD expression" as a denominator — reference share is '
                    'ref / Σ ALL isoforms = "total gene expression"')

    # C: section back-references
    for m in re.finditer(r"§|\bSection\s+\d|as described in section", text, re.I):
        warns.append(f'section back-reference @ …{text[max(0,m.start()-15):m.start()+15]}… '
                     '(legends are standalone)')

    # B: British spellings
    for m in re.finditer(r"\b(colour|coloured|centre|centred|grey|analyse|analysed|"
                         r"normalise|normalised|characterise|characterised|behaviour|"
                         r"labelling|modelling)\b", text, re.I):
        errs.append(f'British spelling "{m.group(0)}" (use US spelling)')

    # C: SF title sentence-case heuristic (WARN)
    tm = re.search(r"\*\*\s*(SF\d+[^|]*)\|\s*([^*]+?)\.?\*\*", text)
    if tm:
        title = tm.group(2).strip()
        words = [w for w in re.split(r"[\s/]+", title) if w]
        sig = [w for w in words if w.lower() not in SMALL and w[:1].isalpha()]
        lc = [w for w in sig if w[0].islower()]
        if sig and len(lc) / len(sig) > 0.4:
            warns.append(f'title may be sentence-case (Title Case required): "{title}" '
                         f'[{len(lc)}/{len(sig)} content words lowercase]')

    return errs, warns


def main():
    base = Path(__file__).resolve().parents[1]
    roots = [Path(p) for p in sys.argv[1:]] or \
            [base / "SupplementalFigures", base / "multipanel"]
    files = sorted({f for r in roots for f in
                    ([r] if r.is_file() else r.rglob("*legend*.md"))})
    if not files:
        print("no *legend*.md found under:", *roots); return 0
    n_err = n_warn = 0
    for f in files:
        errs, warns = check(f)
        if not errs and not warns:
            continue
        rel = f.relative_to(Path.cwd()) if str(f).startswith(str(Path.cwd())) else f
        print(f"\n{rel}")
        for e in errs:
            print(f"  ERROR  {e}"); n_err += 1
        for w in warns:
            print(f"  warn   {w}"); n_warn += 1
    print(f"\n{'='*60}\n{len(files)} legend(s) scanned · {n_err} ERROR · {n_warn} warn")
    if n_err == 0:
        print("mechanical checks clean — now do §A number-matching by hand.")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
