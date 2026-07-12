# SF37 — End-of-window channel decomposition (internal analysis notes)

**Not a manuscript legend.** Internal analysis notes to accompany SF37, following up on
the visible elevation of the summed |SHAP| curve at the far right of the stop-codon
window (Panel B) and, to a lesser extent, at the far left of the AUG window (Panel A).

## Data

Table: `data/sf37_edge_channel_decomposition.tsv` (90 rows).

Per-channel mean |SHAP| aggregated into five position zones within each 500-nt window:

- `left_edge`  — the outermost 50 nt on the left  (positions -250 … -201)
- `left_mid`   — the next 100 nt              (positions -200 … -101)
- `middle`     — the central 200 nt                (positions -100 …  +99)
- `right_mid`  — the next 100 nt                   (positions +100 … +199)
- `right_edge` — the outermost 50 nt on the right  (positions +200 … +249)

Columns: `window` (`atg` / `stop`), `zone`, `channel` (`A`, `C`, `G`, `T`, `frame_0`,
`frame_1`, `frame_2`, `junction`, `rolling_gc`), `nmd_mean_abs_shap`,
`ctrl_mean_abs_shap`, `n_positions`.

Underlying |SHAP| ancestry: pooled across 5 DeepSHAP runs (500 background samples each)
on the atg500_stop500 configuration of the v5_4ct model. Same source as SF37 Panels A/B.

## Where the end-of-window elevation is coming from

### Stop-codon window, right edge (positions +200 … +249) — the strongest end-effect

Top drivers by NMD |SHAP| in this zone:

| Channel     | NMD |SHAP| | Control |SHAP| |
|-------------|-----------|----------------|
| **T**       | 0.00846   | 0.00752        |
| G           | 0.00666   | 0.00528        |
| A           | 0.00589   | 0.00436        |
| C           | 0.00575   | 0.00411        |
| rolling_gc  | 0.00197   | 0.00207        |
| junction    | 0.00161   | 0.00054        |
| frame_0/1/2 | 0.0005–0.0006 | 0.0005–0.0006 |

The uracil (T) channel is the dominant driver here, ~30% larger than the next
nucleotide channel (G). This is consistent with the model attending to
polyadenylation infrastructure in the 3′ terminus of NMD-susceptible transcripts:
the canonical AAUAAA hexamer sits ~10–30 nt upstream of the cleavage/poly-A site,
and a downstream U-rich element (the DSE) sits ~20–40 nt downstream of it. Both
elements are T/U-rich and cluster within ~50 nt of the transcript 3′ end.

The junction-channel |SHAP| is also ~3× higher in NMD than Control at this
right edge — consistent with EJC-proximity signal in the immediate 3′ terminus
of PTC-containing transcripts.

Frame channels contribute almost nothing at this edge (~0.0005 for both classes),
so the elevation is **not** a coding-frame or padding artifact.

### AUG window, left edge (positions -250 … -201) — the smaller left-edge elevation

Top drivers by NMD |SHAP| in this zone:

| Channel     | NMD |SHAP| | Control |SHAP| |
|-------------|-----------|----------------|
| **G**       | 0.00252   | 0.00122        |
| C           | 0.00191   | 0.00096        |
| A           | 0.00153   | 0.00081        |
| T           | 0.00146   | 0.00079        |
| rolling_gc  | 0.00133   | 0.00066        |
| frame_0/1/2 | ~0.0007   | ~0.0004        |
| junction    | ~0.00001  | ~0.00001       |

The G channel leads at the left edge of the AUG window, closely followed by C,
and rolling-GC is also elevated. This is consistent with generally GC-rich
sequence in the far-upstream 5′UTR neighborhood the model is scanning. Frame
channels are the *smallest* contributor here (contrary to a plausible initial
guess) — so the left-edge elevation is not a frame-marker artifact.

A caveat: for transcripts whose 5′UTR is shorter than 250 nt (a substantial
share of the cohort), positions in the far-left zone are zero-padded. The
model's attention to the padded region contributes to the |SHAP| signal but
not necessarily to biological signal; the shape of the class difference in
the padded region is not itself interpretable.

## Bottom line

- The right-edge stop-window elevation is a real signal: T-channel dominance is
  consistent with poly-A signal infrastructure that concentrates at the 3′ end of
  the transcript.
- The left-edge AUG-window elevation is milder, has nucleotide-composition
  origins (G/C leading), and is confounded by zero-padding in short-5′UTR
  transcripts.
- Neither edge is driven by the frame channels.

## Provenance

- Rendered from: `figures/SupplementalFigures/SF37_ShapAcrossWindows/data/sf37_edge_channel_decomposition.tsv`
- Underlying |SHAP| ancestry: `model:09b_export_subgroup_profiles.py` from the 5-run
  DeepSHAP array (`slurm_deepshap_seq_500bg.sh`) at model repo commit `5c19591`.
- Cohort: SF37 held-out test set — NMD n = 2,268, Control n = 7,863.
