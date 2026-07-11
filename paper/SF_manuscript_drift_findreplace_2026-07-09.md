# Manuscript find/replace pairs — SF drift fixes (2026-07-09)

Six drifts flagged by the pass-3 independent reviewer between the Google Doc
manuscript and the current SF24–SF42 legends. All are manuscript-side fixes:
the SFs represent current canonical state (post-TD2-bias fix + SF32/33/34/35
split); the manuscript prose was written before those changes.

Format below is `FIND` (verbatim from the manuscript) → `REPLACE`.
Recommended order: apply top-to-bottom; each change is self-contained.

## Note before starting

The reviewer's Drift #5 ("Fig 4 caption Panel C mislabels the 1,050 group as
`NMD+/PTC−`") was a **false positive**. The current manuscript correctly reads
`1,050 (90%) had a downstream-EJC stop (NMD+/PTC+)`. Skipping.

---

## Drift #3 — Methods: "hidden-PTC subset" → "occult-PTC subset"

Standardize to "occult" (the same paragraph mixes both terms; every SF uses
"occult"; the `feedback_nmd_paper_check_cds_source` memory rule confirms
"occult PTC" as the paper's canonical form).

**FIND:**
```
a hidden-PTC subset
```
**REPLACE:**
```
an occult-PTC subset
```

If any other `hidden-PTC` occurrences exist elsewhere in the doc, apply the
same replacement globally.

---

## Drift #1 — §4 last paragraph: SF35 miscited on the GENCODE-restricted
sentence (should be SF32)

The joint SF32+SF35 figure was split into two separate SFs; this citation was
not updated. SF32 is the GENCODE-restricted cohort at n = 190; SF35 is the
reference AUG-traceable cohort at n = 1,166.

**FIND:**
```
NMD+/PTC+ isoforms and the other groups in the GENCODE restricted set (SF35). In the
expanded set, NMD+/PTC+ isoforms had significantly longer isoforms than both other groups
```
**REPLACE:**
```
NMD+/PTC+ isoforms and the other groups in the GENCODE restricted set (SF32). In the
expanded set, NMD+/PTC+ isoforms had significantly longer isoforms than both other groups (SF35)
```

---

## Drift #4 — §4 paragraph 4: 99% (487/492) TD2-downstream claim is uncited

The 99% (487/492) and 78% (384/492) claims both come from SF34 (occult-PTC
subset). Only the second is currently cited. Add SF34 to the first.

**FIND:**
```
99% (487/492) of the TD2-called CDS were
```
**REPLACE:**
```
99% (487/492) of the TD2-called CDS were downstream from the reference AUG (SF34)
```

If the above wraps awkwardly, the alternative one-shot form is to move both
into a single citation later in the sentence — but the reviewer's preferred
minimum is the added `(SF34)` after the "downstream from the reference AUG"
phrase.

---

## Drift #2 + Drift #6 — §4 paragraph 3: SF30/SF31 sit inside an n = 190
paragraph but are now n = 1,166 / n = 1,050 (ref-AUG anchor)

SF30 and SF31 were re-anchored from TD2 stop → reference AUG-traced stop
this session, and the population dropped from n = 2,790 / 1,211 → n = 1,166 /
1,050. The n = 190 GENCODE-restricted paragraph continues to cite them
without noting the scope shift. Fix in place with a one-clause insertion.

**FIND:**
```
We observed a strong enrichment of stop codons far
upstream of the terminal exon junction in NMD isoforms (Fig 3D - PTCs in NMD) with a clear
dose-response relationship between upstream distance of the stop codon and the magnitude of
the NMD response (SF30 - PTC distance dose response). There was a similar relationship
with the number of downstream EJCs, where the magnitude of the NMD response peaked at
4-5 downstream EJCs from the PTC (SF31 - NMD Effect Size by Number of EJCs).
```
**REPLACE:**
```
We observed a strong enrichment of stop codons far
upstream of the terminal exon junction in NMD isoforms (Fig 3D - PTCs in NMD). Extending
the same analysis to the broader reference AUG-traceable cohort (n = 1,166), we saw a clear
dose-response relationship between the distance of the reference-anchored stop codon to
the last EJC and the magnitude of the NMD response (SF30 - PTC distance dose response).
There was a similar relationship with the number of downstream EJCs, where the magnitude
of the NMD response peaked at 4–5 downstream EJCs from the PTC (SF31 - NMD Effect Size
by Number of EJCs).
```

This single edit addresses both Drift #2 and Drift #6 (SF30 anchor/scope
signaling missing from Fig 3D framing).

---

## Optional but recommended: numerical drift already noted elsewhere

`results_to_code_map.md` at line ~626 still describes SF30 as `n = 2,790` — a
pre-fix number. If someone reads that map alongside the manuscript, they will
be confused. Repo-side change; not manuscript.

---

## Summary of citations after the pairs above are applied

| Sentence | Was | Now |
|---|---|---|
| §4 GENCODE-restricted 3′UTR result | (SF35) | (SF32) |
| §4 expanded-set 3′UTR result | (unlabeled follow-on) | (SF35) |
| §4 SF30/SF31 introduction | in n = 190 paragraph, no scope caveat | in a scope-tagged sentence anchored to n = 1,166 |
| §4 99% TD2-downstream | uncited | (SF34) |
| Methods "hidden-PTC subset" | hidden-PTC | occult-PTC |

Five sentences touched. No changes to numbers, panel counts, or any other
prose.
