# Deposit provenance — INTERNAL RECORD

> ⚠️ **Internal only. Do not copy into paper materials.**
> Not in `nmd_lung_longread_2026`, not in the Zenodo record, not in the manuscript or
> supplement. Pete, 2026-07-25: *"We should record this for ourselves, but not in the paper
> materials."*

## Why this file exists

The project's data **is** the four cell types — AT2, LAE, FB, MV. That is what the paper
reports, what GEO GSE329233 carries, and what the Zenodo source-data record contains. None of
those materials should describe the dataset as a subset of anything larger, because from the
paper's point of view it isn't.

Internally, though, we should be able to answer "how exactly was the deposited count matrix
produced, and is it faithful?" without re-deriving it. That is what this note is for.

## What was done

Sequencing originally covered additional culture conditions beyond the four cell types in this
paper. The deposited count matrices were assembled to contain **only** the 26 samples of the
four studied cell types — 13 donor-matched DMSO / SMG1i pairs — with columns named
`<donor>_<treatment>_<cell type>` using canonical names (AT2, LAE, FB, MV).

Selecting at source, rather than loading everything and filtering downstream, is what allows
the analysis code to be written natively for four cell types (REPRODUCIBILITY_PLAN §2b,
"Option C"). No deposited or shipped script drops samples.

## The scripts

`deposit/build_4ct_count_matrix.R` and `deposit/build_4ct_salmon_counts.R` (this repo,
internal). They are **deliberately not shipped** — they describe how a data package was
assembled, not how any published result was produced, and their inputs are not deposited, so a
reader could not run them anyway.

Non-obvious detail worth keeping: the short-read salmon matrix could **not** be subset by
column-name prefix. Three of the `DD`-prefixed columns are in fact ALI cultures (donors 001V,
027U, 029T), the short-read cell-type code is `AT2` where the long-read code is `AT`, and two FB
columns carry a `_clean` suffix. A naive prefix filter silently retains six wrong samples. The
script therefore takes its sample list from the published SR DGEList itself.

## Faithfulness

Both scripts self-verify and abort on mismatch. Confirmed 2026-07-24:

| output | check | result |
|---|---|---|
| `nmd_lungcells_counts_4ct.csv` | counts vs `dge_isoform_longread_2026.3.3` | **identical** |
| `salmon_gene_counts_4ct.csv` | counts vs `dge_shortread_gene_2026.3.2` | **identical** |

**The deposit is trimmed to observed features.** Isoforms 645,272 → 614,992 and genes
78,899 → 46,571; a feature with zero counts in every sample is not part of this dataset, and
leaving such rows in invites questions the paper does not answer.

Verified result-neutral before trimming: the `filterByExpr` sets are identical from the full and
trimmed universes, and percent-output-lost is unchanged to four decimals. The published filtered
sets contain no all-zero features.

Trimmed consistently across every keyed file — both count matrices, SQANTI classification,
corrected FASTA, filtered GTF and corrected CDS GFF3 — all now carrying an identical 614,992
isoform ID set. `trim_deposit_to_observed.sh` (this repo, internal) performs it; the GTF/GFF3
pass is two-pass so original line order is preserved and gene records are retained wherever a
gene keeps at least one transcript.

*(Superseded: an earlier version of this note argued against trimming, on the grounds that
trimming the counts alone would leave the deposit's files disagreeing. That was resolved by
trimming all of them.)*

## Open action — GEO

**GEO GSE329233 currently includes the other culture conditions and they must be removed**, so
the record matches the paper and the Zenodo deposit. Pete, 2026-07-25: do this **after
verification is complete** — the verification work runs against local data and does not depend
on the GEO record, so removing samples earlier would gain nothing and risks disturbing a
submission mid-flight. Submission files: `~/claude_projects/ncbi_submissions/nmd_lung_cells/`.

## Two isoform count matrices (2026-07-25)

The published analyses did **not** all start from the same isoform universe:

| branch | matrix | isoforms |
|---|---|---|
| §1/§2 gene- and isoform-level mashr | SQANTI3-filtered | 645,272 |
| §3/§4 Isopair | isocall, unfiltered by SQANTI | 645,273 |

Evidence: `dge_isoform_longread_2026.3.3.rds` has 645,272 rows and lacks `ENST00000441168.6`;
`expression_data.rds` (the post-filter Isopair CPM matrix) **contains** it.

`ENST00000441168.6` (*ANKMY1*, `ENSG00000144504.17`) carries a 1-bp exon, which is why SQANTI
drops it. It is not incidental to §3/§4: it clears both the 5% isoform-proportion filter and
`filterByExpr`, and sits in the published 95,623-isoform Isopair analysis universe.

Depositing only the filtered matrix covered 95,622 / 95,623 of that universe. Substituting it
would not have been inert — removing the isoform also removes it from its gene's
isoform-proportion denominator, moving sibling `ENSG00000144504.17.novel34` across the
`>= 0.05` boundary. A two-way (−1, +1) change to the universe, not a one-sided drop.

Pete, 2026-07-25: deposit both (option **b**), so each branch reads the matrix it actually
used. `build_4ct_isocall_counts.R` (this repo, internal) builds the second one and verifies:
every published Isopair isoform present (0 missing), `ENST00000441168.6` present and non-zero
(21 counts across the 26 samples), and counts for all 614,992 shared isoforms **identical** to
the filtered matrix.

Both are trimmed to features observed in the 26 samples — 614,992 and 614,993 rows, differing
by exactly that one isoform.

*(Superseded: an earlier section here argued the deposit files should not be trimmed to
observed isoforms. They were trimmed; the section contradicted §Faithfulness above and has
been removed.)*
