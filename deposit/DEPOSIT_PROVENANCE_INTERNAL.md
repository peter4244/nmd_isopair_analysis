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

All 645,272 isoforms are retained in the isoform matrix, including those with zero counts in
these samples — the percent-transcriptional-output-lost measure and the
cell-type-restricted-expression analyses are defined over the full universe, so trimming would
change published results.

## Consistency check still open

**GEO GSE329233 must contain only the four cell types' reads.** If the submitted record includes
other conditions, it will not match the paper or the Zenodo deposit, and a reader comparing them
would see a discrepancy. Verify against
`~/claude_projects/ncbi_submissions/nmd_lung_cells/` before submission.
