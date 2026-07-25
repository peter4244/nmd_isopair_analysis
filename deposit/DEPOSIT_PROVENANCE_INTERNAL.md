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

All 645,272 isoforms are retained in the isoform matrix. (An earlier version of this note
claimed trimming would change published results. That was wrong — it was tested and is lossless.
The actual reasons for retaining them are below.)

## Open action — GEO

**GEO GSE329233 currently includes the other culture conditions and they must be removed**, so
the record matches the paper and the Zenodo deposit. Pete, 2026-07-25: do this **after
verification is complete** — the verification work runs against local data and does not depend
on the GEO record, so removing samples earlier would gain nothing and risks disturbing a
submission mid-flight. Submission files: `~/claude_projects/ncbi_submissions/nmd_lung_cells/`.

## Why the deposited files were not trimmed to observed isoforms

Considered and rejected 2026-07-25. The isoform matrix carries 645,272 isoforms, of which 30,280
have no counts in these samples. Trimming to the 614,992 observed was **verified lossless** —
identical filterByExpr set, and percent-output-lost identical to four decimals — but rejected
because:

1. The classification, GTF, GFF3 and FASTA all carry the same 645,272 isoforms. Trimming the
   counts alone would leave the deposit's files disagreeing with each other, which is a worse
   inconsistency than the one it fixes. Trimming all of them means parsing 3.8 GB of
   GTF/GFF3/FASTA for no analytical gain.
2. SQANTI filters on structural criteria rather than expression, so a filtered set containing
   isoforms with no assigned counts is ordinary and not, by itself, informative about anything.

The README note that drew attention to it ("including those with zero counts") was removed
instead. Had we trimmed, exactly two baseline lines would have changed — the DGEList and CPM
dimension printouts in `verify_section3_p1.txt` — with all nine numeric-result lines unaffected.
