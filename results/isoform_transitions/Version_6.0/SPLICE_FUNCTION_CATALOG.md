# Event Detection Function Catalog

Reference for all functions in `scripts/event_detection_functions.R`, organized by
the aspect of splicing biology they address.

## Event Type Reference

| Event Type | Splicing Mechanism | Emitting Function(s) | Detection Step |
|---|---|---|---|
| **Alt_TSS** | Alternative promoter / first exon | `detect_events_for_pair` Step 3, via `detect_tss_change` | Terminal |
| **Alt_TES** | Alternative polyadenylation / last exon | `detect_events_for_pair` Step 3, via `detect_tes_change` | Terminal |
| **A5SS** | Alternative 5' splice site (donor) | `detect_events_for_pair` Step 2b, via `detect_shared_boundary_event` / `check_boundary_within_exon` | Boundary |
| **A3SS** | Alternative 3' splice site (acceptor) | `detect_events_for_pair` Step 2b, via `detect_shared_boundary_event` / `check_boundary_within_exon` | Boundary |
| **Partial_IR_5** | Partial intron retention at donor | `detect_events_for_pair` Step 2b, via `detect_shared_boundary_event` | Boundary |
| **Partial_IR_3** | Partial intron retention at acceptor | `detect_events_for_pair` Step 2b, via `detect_shared_boundary_event` | Boundary |
| **IR** | Full intron retention | `detect_events_for_pair` Step 2a, via `detect_ir_simple` | IR |
| **IR_diff_5** | Intron retention + 5' boundary mismatch | `detect_events_for_pair` Step 2a, via `detect_ir_simple` | IR |
| **IR_diff_3** | Intron retention + 3' boundary mismatch | `detect_events_for_pair` Step 2a, via `detect_ir_simple` | IR |
| **IR_diff_5_3** | Intron retention + both boundaries mismatch | `detect_events_for_pair` Step 2a, via `detect_ir_simple` | IR |
| **SE** | Exon skipping (cassette exon) | `detect_events_for_pair` Step 2c (inline flanking check) | Skipping |
| **Missing_Internal** | Internal exon absent (flanking not met) | `detect_events_for_pair` Step 2c + non-overlap gap zone | Skipping |

All events carry a **direction** (GAIN or LOSS) from the comparator's perspective.

---

## Function Catalog

### Transcript orientation

**`order_exons_biological(exons, strand)`** — Sorts exons in transcription order
(5' to 3'). This establishes the biological frame of reference: RNA polymerase reads
the template strand from 5' to 3', so on the plus strand that's low-to-high
coordinates, and on minus strand it's high-to-low. Every downstream function depends
on this ordering to correctly identify which end is TSS vs TES and which splice site
is donor vs acceptor. *Does not emit events.*

---

### Transcription initiation and termination

These detect differences at the two ends of the transcript that are **not determined
by the spliceosome** — they reflect alternative promoter usage (TSS) or alternative
polyadenylation (TES).

**`detect_tss_change(exon_dom, exon_non_dom, strand)`** — Compares the 5' boundary
of each isoform's first exon. A difference beyond 20bp indicates the two isoforms
originate from different promoters (alternative first exon). The tolerance accounts
for imprecise 5' mapping in long-read data. Additionally, even if the TSS coordinate
difference is within tolerance, returns TRUE when the first exons don't overlap at
all — this catches cases where a short extra terminal exon creates a small coordinate
difference but represents a genuine structural change. *Gate for* -> **Alt_TSS**.

**`detect_tes_change(exon_dom, exon_non_dom, strand)`** — Same logic at the 3' end.
A difference beyond 20bp indicates alternative polyadenylation or an alternative last
exon. Also returns TRUE when last exons don't overlap, regardless of coordinate
distance. *Gate for* -> **Alt_TES**.

**`compute_missing_terminal_exons_tss(dom_exons, comp_exons, strand)`** — When one
isoform's promoter is further upstream, this walks from that promoter inward,
collecting exonic regions (skipping introns) that the shorter isoform lacks. These
are the exons that exist only because a different promoter was used — the spliceosome
never sees them in the shorter isoform because transcription didn't extend that far.
*Populates `missing_terminal_exons` field on* -> **Alt_TSS**.

**`compute_missing_terminal_exons_tes(dom_exons, comp_exons, strand)`** — Mirror
image at the 3' end. Collects exonic regions downstream of the shorter isoform's
polyadenylation site. *Populates `missing_terminal_exons` field on* -> **Alt_TES**.

**`compute_orphan_terminal_exons_tss(walk_ordered, ref_ordered)`** — Finds exons at
the 5' end of one isoform that sit within the other isoform's genomic span but don't
overlap any of its exons. These represent a distinct first exon choice — the
transcript "lands" in a completely different genomic region after the promoter,
suggesting an alternative first exon event rather than just a shifted TSS.
*Populates `orphan_terminal_exons` field on* -> **Alt_TSS**.

**`compute_orphan_terminal_exons_tes(walk_ordered, ref_ordered)`** — Same logic at
the 3' end. An alternative last exon that doesn't overlap the other isoform's
terminal exons. *Populates `orphan_terminal_exons` field on* -> **Alt_TES**.

---

### Splice site choice (spliceosome-mediated)

These detect differences in where the spliceosome cuts — the core of alternative
splicing.

**`detect_shared_boundary_event(exon_dom, exon_non_dom, strand, ...)`** — The main
splice-site classifier. Compares two overlapping exons from different isoforms by
examining their donor (5' splice site) and acceptor (3' splice site) boundaries.
*Emits (via return to Step 2b):*

- **One boundary shared, one differs**: The spliceosome used the same site on one end
  but chose a different site on the other. Classified by the size of the shift:
  - <100bp -> **A5SS** (alternative donor) or **A3SS** (alternative acceptor). These
    represent nearby competing splice sites where the spliceosome "chose" one over the
    other.
  - >=100bp -> **Partial_IR_5** or **Partial_IR_3** (partial intron retention). The
    extension is large enough that it likely represents retention of part of an intron
    rather than selection of a nearby competing splice site.

- **Both boundaries differ (internal exons)**: Two independent splice-site choices
  happened at the same exon. Decomposed into primary event + `second_event` (one per
  boundary), each independently -> **A5SS**/**Partial_IR_5** +
  **A3SS**/**Partial_IR_3**.

- **Both boundaries differ (terminal exon involved)**: Only the internal-facing
  boundary is emitted as a splice-site event. The outer boundary is handled by
  Alt_TSS/Alt_TES. When asymmetric (one exon is terminal, other is not), the
  splice-site-facing boundary is emitted as `second_event`.

**`check_boundary_within_exon(comp_boundary, dom_boundary, ...)`** — Fallback
detection for exon pairs that overlap but don't share an exact boundary. Determines
whether the comparator's splice site falls within the dominant's exon (the dominant
exon is longer -> LOSS) or extends beyond it (the comparator exon is longer -> GAIN).
Uses the same <100bp vs >=100bp threshold. *Returns classification to
`detect_shared_boundary_event`* -> **A5SS**, **A3SS**, **Partial_IR_5**,
**Partial_IR_3**, or **IR** (when extension spans into flanking exon).

**`check_spans_flanking_exons(exon_dom, exon_non_dom, ...)`** — Distinguishes partial
intron retention from full intron retention. If the longer exon's extension reaches
into a neighboring exon in the other isoform, it's spanning an entire intron (IR),
not just partially retaining one. This is the biological distinction between "the
spliceosome failed to remove this intron" (IR) vs "the spliceosome used a slightly
different splice site" (Partial_IR). *Returns boolean to
`detect_shared_boundary_event`; suppresses Partial_IR emission when TRUE (defers to
IR detection in Step 2a).*

---

### Intron retention

**`detect_ir_simple(exon, other_exons)`** — Tests whether a single exon in one
isoform spans two or more exons in the other. This is the hallmark of intron
retention: one isoform's spliceosome excised the intron (producing two separate
exons), while the other isoform retained it (one continuous exon spanning the
intronic region). *Gate for Step 2a* -> **IR**, **IR_diff_5**, **IR_diff_3**,
**IR_diff_5_3** (subtype classification is inline in `detect_events_for_pair`).

The IR_diff subtypes capture cases where the retained-intron exon's boundaries don't
exactly match the outermost split exons — indicating the retention event co-occurs
with a boundary shift at the 5' end, 3' end, or both.

---

### Exon skipping

**`detect_se(comparison)`** — (Standalone/legacy version, not called by
`detect_events_for_pair`.) Identifies cassette exons from a union-exon comparison
dataframe. *Would emit* -> **SE**. In the current pipeline, SE detection is performed
inline in Step 2c of `detect_events_for_pair`.

**Step 2c inline logic** (within `detect_events_for_pair`) — For each exon present in
only one isoform (within the other's genomic span, no overlap), checks strict
flanking: both neighboring exons must overlap the other isoform. *Emits:*

- **SE** — flanking condition met. Classic cassette exon: the spliceosome either
  includes or excludes an entire exon, and the flanking exons confirm the shared
  pre-mRNA context.
- **Missing_Internal** — flanking condition not met. The exon is absent but the
  structural context is more complex (e.g., terminal exon of one isoform, or
  neighbors also differ). Also emitted for non-overlapping isoforms' gap-zone dom
  exons.

---

### Splice junction tracking

These don't detect events but record the splice junctions involved, providing
evidence for each event.

**`compute_junctions(exons_ordered)`** — Derives splice junctions from consecutive
exon pairs. Each junction represents one intron removal by the spliceosome (donor
site of exon N -> acceptor site of exon N+1). *Does not emit events.*

**`junctions_touching_range(jxn_vec, range_start, range_end)`** — Finds junctions
where at least one splice site falls within the event range. Used for A5SS/A3SS/SE
events where the affected splice site is at the event boundary. *Populates
`dom_junctions`/`comp_junctions` fields on* -> **A5SS**, **A3SS**, **Partial_IR_5**,
**Partial_IR_3**, **SE**, **Missing_Internal**.

**`junctions_within_range(jxn_vec, range_start, range_end)`** — Finds junctions fully
contained within a range. Used for IR events: the dominant's split exons have
junctions inside the retained-intron region. *Populates
`dom_junctions`/`comp_junctions` fields on* -> **IR**, **IR_diff_\***.

**`junctions_spanning_range(jxn_vec, range_start, range_end)`** — Finds junctions
that skip over a range entirely. Used for SE events: the skipping isoform has a
single junction that leaps over the cassette exon. *Populates
`dom_junctions`/`comp_junctions` fields on* -> **SE**, **Missing_Internal**.

**`format_junctions(jxn_vec)`** — Utility to serialize junctions as a comma-separated
string. *Does not emit events.*

---

### Isoform-level guards

**`check_isoform_overlap(comp_ordered, dom_ordered)`** — Tests whether the two
isoforms share any exonic overlap at all. If not, the internal splicing comparison
(IR, A5SS, A3SS, SE) is skipped — there's no shared pre-mRNA context where the
spliceosome made different choices. Only terminal events (Alt_TSS/Alt_TES) and
gap-zone Missing_Internal are meaningful for non-overlapping isoforms. *Does not emit
events; controls flow into Step 2 vs gap-zone path.*

---

### Orchestrator

**`detect_events_for_pair(dominant_exons, comparator_exons, ...)`** — Runs the full
detection pipeline for one isoform pair. All events are ultimately emitted here.
Detection proceeds in biological priority order:

| Step | What it detects | Event types emitted |
|---|---|---|
| 2a | Intron retention | IR, IR_diff_5, IR_diff_3, IR_diff_5_3 |
| 2b | Splice-site boundary shifts | A5SS, A3SS, Partial_IR_5, Partial_IR_3 |
| 2c | Exon skipping / missing exons | SE, Missing_Internal |
| gap | Non-overlapping dom exons in comp span | Missing_Internal |
| 3 | Terminal boundary differences | Alt_TSS, Alt_TES |
