# GEO SuperSeries submission package — NMD lung cells

This directory contains everything needed to build the GEO SuperSeries submission
for the NMD lung-cell study (short-read + long-read RNA-seq). Most of the
package is metadata + scripts that run on the cluster to stage the actual files
from their production locations.

## Package layout

```
geo_submission/
├── README.md                              (this file)
├── sample_metadata_longread.csv           38-sample phenotype table (long-read)
├── metadata/
│   ├── series.md                          SERIES text (title, abstract, overall design)
│   ├── protocols.md                       PROTOCOLS + DATA PROCESSING text
│   ├── samples_longread.tsv               SAMPLES section rows (long-read, 38 rows)
│   └── samples_shortread.tsv.template     SAMPLES template — filled on cluster
└── scripts/
    ├── stage_longread.sh                  cluster: stage uBAMs + isocall outputs
    └── stage_shortread.sh                 cluster: stage FASTQs + salmon matrix
```

After the staging scripts run on the cluster, a parallel tree is created under
`/proj/…/geo_submission/`:

```
geo_submission/                    (on cluster — produced by the staging scripts)
├── longread/
│   ├── raw/                       {alias}.bam symlinks → bams/ALL/
│   ├── processed/                 isocall_transcript_counts.tsv, isocall_isoforms.gtf.gz, sample_metadata_longread.csv
│   ├── md5sums_raw.txt
│   ├── md5sums_processed.txt
│   └── metadata/
│       ├── raw_files_longread.tsv          paste into RAW FILES section
│       └── processed_files_longread.tsv    paste into PROCESSED DATA FILES section
└── shortread/
    ├── raw/                       {alias}_R1.fastq.gz / {alias}_R2.fastq.gz symlinks
    ├── processed/                 salmon_gene_counts_length_scaled.tsv, sample_metadata_shortread.csv
    ├── md5sums_raw.txt
    ├── md5sums_processed.txt
    └── metadata/
        ├── samples_shortread.tsv           paste into SAMPLES section
        ├── raw_files_shortread.tsv         paste into RAW FILES section
        ├── processed_files_shortread.tsv   paste into PROCESSED DATA FILES section
        └── paired_end_shortread.tsv        paste into PAIRED-END section
```

## Submission structure

One **SuperSeries** with two **SubSeries**:

| SubSeries  | Platform                         | Raw files           | Processed files |
|------------|----------------------------------|---------------------|-----------------|
| Short-read | Illumina paired-end              | FASTQ (R1+R2)       | `salmon_gene_counts_length_scaled.tsv` (nf-core/rnaseq, STAR + Salmon v1.10) |
| Long-read  | PacBio HiFi Iso-Seq (Kinnex)     | uBAM, one per sample | `isocall_transcript_counts.tsv`, `isocall_isoforms.gtf.gz` (isocall v0.15.0) |

## End-to-end workflow

### 1. Fill in the metadata spreadsheet (local)

1. Download the GEO high-throughput sequencing template from NCBI:
   <https://www.ncbi.nlm.nih.gov/geo/info/seq.html> → `seq_template.xlsx`.
2. Make three copies:
   - `seq_template_superseries.xlsx`
   - `seq_template_shortread.xlsx`
   - `seq_template_longread.xlsx`
3. For each, paste the relevant text blocks from `metadata/series.md` and
   `metadata/protocols.md` into the SERIES and PROTOCOLS sheets.
4. Leave the SAMPLES, RAW FILES, PROCESSED DATA FILES, and PAIRED-END sections
   empty for now — the cluster scripts in step 2 will produce TSV blocks that
   paste directly into those sections.
5. Resolve every `TODO` placeholder in `metadata/series.md` and
   `metadata/protocols.md` (treatment dose/duration, RNA kit, library kit,
   instrument model, software versions). *(isocall default parameter values
   are already pinned from `isocall call --help` output.)*

### 2. Stage the files (cluster)

SSH to the cluster, then:

```bash
cd <repo checkout>/geo_submission/scripts

bash stage_longread.sh       # symlinks uBAMs, copies isocall outputs, md5s
bash stage_shortread.sh      # symlinks FASTQs, merges & copies Salmon matrix, md5s
```

Each script is idempotent — re-runnable without side effects. If any input
file is missing, the script exits with a clear error listing which samples
or paths are unresolved; fix and rerun.

Before running, edit these lines in each script if your paths differ:

- `stage_longread.sh`: `PROJECT_REPO`, `UBAM_SOURCE_DIR`, `STAGING_ROOT`,
  `ISOCALL_COUNT_MATRIX`, `ISOCALL_ISOFORMS_GTF`.
- `stage_shortread.sh`: `PROJECT_REPO`, `PHENO_FB_MV`, `PHENO_DD_DO_AT`,
  `FASTQ_DIRS`, `COUNTS_BATCH1`, `COUNTS_BATCH2`, `INSTRUMENT_MODEL`,
  `READ_LENGTH`.

### 3. Paste TSV blocks into the spreadsheets (local)

Pull the four TSVs (SAMPLES, RAW FILES, PROCESSED DATA FILES, PAIRED-END) out
of each SubSeries' `metadata/` directory on the cluster and paste them into the
corresponding sections of the short-read and long-read seq_template.xlsx files.

### 4. Upload to GEO

Request an Aspera/FTP upload slot from GEO (<geo@ncbi.nlm.nih.gov>) using your
existing GEO account. You will get a unique host + token. Upload both
SubSeries file trees under the same submission:

```bash
# Example pattern — exact command comes from GEO email
ascp -QT -l 300m -k 1 -d \
  longread/raw/  longread/processed/ \
  shortread/raw/ shortread/processed/ \
  <user>@<geo-host>:/<path>
```

Include **all** raw files, all processed files, both MD5 files, and both
completed XLSX spreadsheets in the upload.

### 5. Submit the SuperSeries via the GEO web form

- Log into GEO with the existing account.
- Create a new SuperSeries submission; link the two SubSeries GSE numbers once
  they are assigned by GEO.
- Paste the SuperSeries-level SERIES text from `metadata/series.md`.

## Sample aliases

Every sample in this submission uses the alias pattern
`{cell_type}_{donor}_{treatment}_{SR|LR}` — e.g.,
`DD_ALI_027U_Smg1i_LR`, `MV_001V_DMSO_SR`. The alias is the primary key linking:

- `sample_metadata_longread.csv` / `sample_metadata_shortread.csv`
- The SAMPLES row in the GEO spreadsheet
- The staged raw filename (`{alias}.bam`, `{alias}_R1.fastq.gz`, `{alias}_R2.fastq.gz`)

Cell-type codes used in the alias:

| Code    | Description |
|---------|-------------|
| `AT`    | Primary alveolar type 2 cells (AT2) |
| `DD`    | Differentiated airway epithelial cells, submerged |
| `DD_ALI`| Differentiated airway epithelial cells, air-liquid interface |
| `DO`    | Undifferentiated donor airway epithelial cells |
| `FB`    | Lung fibroblasts |
| `MV`    | Lung microvascular endothelial cells |

## Before you submit — checklist

- [ ] All `TODO_` placeholders in `metadata/series.md` and `metadata/protocols.md` resolved.
- [ ] Software versions confirmed (nf-core/rnaseq, STAR, Salmon, minimap2, isocall).
- [ ] Illumina instrument model confirmed (e.g., NovaSeq X / NovaSeq 6000 / NextSeq).
- [ ] PacBio instrument confirmed (Revio — inferred from `R84050` run ID prefix).
- [ ] Kinnex / MAS-Seq library prep confirmed with sequencing core.
- [ ] SMG1i dose + treatment duration in the Overall design paragraph.
- [ ] Cell-culture medium details per cell type in the Growth protocol.
- [ ] Contributor list finalized (currently: Castaldi, Peter, J).
- [ ] `stage_longread.sh` reports no MISSING samples.
- [ ] `stage_shortread.sh` reports no MISSING samples.
- [ ] Both MD5 files generated and match uploaded files.
- [ ] Processed count matrix columns match sample aliases in `sample_metadata_*.csv`.

## Notes and caveats

- The long-read uBAMs in the pheno file reference per-batch raw directories
  (`.../data/raw/{20250804,20250909,20251112}/…/IsoSeqX_bc*.bam`). The
  staging script assumes all 38 uBAMs have been consolidated under
  `.../data/raw/bams/ALL/` with the same basenames. If the basenames differ
  in the consolidated directory, update the `source_bam_basename` column of
  `sample_metadata_longread.csv` or edit `stage_longread.sh` accordingly.

- Short-read FASTQ discovery in `stage_shortread.sh` uses a small list of
  common Illumina naming patterns. If the sequencing core used a different
  pattern, the script will flag each unresolved sample and exit — adjust the
  `find_fastq()` helper as needed.

- Two batches of short-read runs exist (`20250811` and `20260129`). The
  staging script searches both. The two corresponding Salmon count matrices
  (`20250718` and `20250811` under `nfcore-rnaseq/`) are merged into a
  single matrix on disk (outer join on `gene_id`) before staging.

- GEO encourages raw data in its original unaligned form for long-read PacBio
  (HiFi uBAMs preserve per-base quality scores and kinetic tags that are lost
  when converting to FASTQ). This submission uses uBAMs.
