# Scientific-Figure Validator: Architecture Plan (v2)

Supersedes `validator_plan_v1.md`. v1's recommendation ("Architecture A": typed
element enumerator + pairwise rule registry + `RenderContext` attached to the
figure, 6-phase / ~7-day migration) is **kept in outline but corrected in four
technical specifics** and **re-sequenced** so the earliest work closes the
currently-open figure defects. The scale-based objection to v1 is rejected (see
§2.0).

All file:line references are to the canonical copy at
`/Users/petecastaldi/claude_projects/nmd/figures/lib/validate_figure_layout.py`
and `.../lib/ggplot_style.py` as of 2026-07-10.

---

## 1. Executive summary

Keep the direction of v1 — a typed enumerator that walks the figure once, a rule
layer that expresses relationships declaratively, and a readability model that
survives composites — but build it on the corrected substrate below:

1. **Rule model is UNARY + PAIRWISE + N-ARY (group) from day one**, not
   pairwise-only. Two of the existing checks are irreducibly n-ary:
   `assert_style_symmetric` (validate:1040) compares a font-size role across *all*
   sibling axes, and the horizontal-alignment check (validate:592, Check 5) groups
   3+ texts by y-band. A `(type_a, type_b) -> rule` table cannot host either
   without faking it.
2. **`Style` is the top-level concept, passed explicitly** as a parameter, never
   attached via `WeakKeyDictionary` or `fig._nmd_render_context`. A `Style` is the
   complete spec a figure conforms to — its *render envelope* (the paper-vs-grant
   scaling constants, formerly `RenderContext`) plus its *visual identity*
   (palette, typography, spine/panel-letter convention). Two shipped instances,
   `PAPER` and `GRANT`, with the type open for future media (POSTER, SLIDE,
   other-journal). The code already threads `native_width_in` explicitly through
   `render_and_validate` → `assert_docx_readable` (ggplot_style:170, validate:1117);
   we generalize that into an explicit `style=` argument, we don't reverse it into
   fig-attached state.
3. **Migration oracle is identical VALIDATOR VERDICTS + metadata-pinned PNGs**,
   not byte-identical PDFs. matplotlib stamps `CreationDate`/`ModDate` into PDF
   metadata (`fig.savefig` at ggplot_style:204–205 passes no `metadata=`), so a
   SHA-256-of-PDF oracle false-diffs on the first re-render.
4. **Coverage-gap detection is a first-class check.** v1's `UnknownElement`
   catches unknown artist *types*; it does not catch a known element *pair or
   group that no rule covers*. We add an explicit coverage assertion so silent
   holes in the registry are visible.

Plus two decisions v1 deferred or omitted:

- **Packaging/SSOT (§4).** A second, **already-drifted** copy of the toolchain
  exists under `~/.claude/utils/`. Recommendation: extract `figures/lib/` into an
  installable, dual-pushed package and delete the stale copy.
- **Tiering that merges with near-term figure work (§5).** Phase 1 is a set of
  **architecture-agnostic in-place fixes** that are correct under any future
  design and that close the open defects (SF31, SF42, the legend-as-data mask,
  the B3 residuals, the SF29/33/34 tick-crowding gaps). The substrate work
  (enumerator, rule model, context, packaging) follows and is pure refactor —
  verdict-preserving by construction.

---

## 2. What changed from v1, and why

### 2.0 Scale reframing (the adversarial critique's scale bucket is REJECTED)

The critique assumed the library serves "15 figures / 2 people" and concluded the
architecture investment is disproportionate. That framing is wrong. This is
**durable cross-project infrastructure**: it continuously services the NMD paper's
supplemental figures (17 gated SFs today — SF25–SF42), the main-text figures, and
the NIH grant figures (`grant_body_fs`/`grant_effective_pt`, ggplot_style:322–340),
and it will service the next paper and the next grant over a multi-year horizon.
The amortization denominator is "every figure this lab ships," not "the 17 SFs in
front of us." The long horizon *justifies* investing in the real substrate; we do
**not** defer the architecture. What we do defer is nothing — we sequence so the
first phase also pays down concrete defects.

### 2.1 N-ary rule model, not pairwise-only (ACCEPTED)

v1's `RULE_REGISTRY: dict[frozenset[type], list[Rule]]` (v1 §2.A) is pairwise. But:

- `assert_style_symmetric` (validate:1040–1100) builds `per_axis[i][role] = size`
  over **all** `fig.axes`, then for each role compares `max(sizes) - min(sizes)`
  against `tolerance_pt`. The unit of judgment is *the whole set of sibling axes
  for one role*, not a pair. Forcing this into pairwise rules means emitting
  O(n²) axis-pairs and losing the natural "3.0pt on axes [1] vs 5.0pt on axes
  [0,2]" report (validate:1091).
- Check 5, horizontal alignment (validate:592–615), sorts texts by `y`, greedily
  partitions into y-bands within `alignment_y_tol`, and for each band of `len ≥ 3`
  reports deviation from the band mean. Also irreducibly a group operation.

**Design.** Three rule kinds, one registry:

```python
# figures/lib/rules/model.py
Severity = Literal["ERROR", "WARNING", "INFO"]

@dataclass(frozen=True)
class UnaryRule:
    applies_to: type[FigureElement]
    fn: Callable[[FigureElement, Style], Issue | None]

@dataclass(frozen=True)
class PairwiseRule:
    a: type[FigureElement]
    b: type[FigureElement]
    fn: Callable[[FigureElement, FigureElement, Style], Issue | None]
    symmetric: bool = True            # if False, order (a,b) is meaningful

@dataclass(frozen=True)
class GroupRule:
    # selector partitions the element list into 0..N groups; fn judges each group.
    selector: Callable[[list[FigureElement]], list[list[FigureElement]]]
    fn: Callable[[list[FigureElement], Style], list[Issue]]
    name: str
```

The runner dispatches by kind (§3.4). `style_symmetric` and `alignment` are
`GroupRule`s; `legend_clearance`, `tick_disjoint`, `text_text_overlap` are
`PairwiseRule`s; `within_canvas`, `docx_readable` are `UnaryRule`s. Worked
conversions for one group rule and one pairwise rule are in §3.5.

### 2.2 `Style` as the top-level concept; explicit, not fig-attached (ACCEPTED + Pete's reframe)

The critique's decision to pass context explicitly stands. Pete's reframe promotes
the abstraction: the top-level concept is a **`Style`** — the complete spec a
figure conforms to — and the paper-vs-grant scaling numbers become the *render
envelope* field inside it. v1 §3 offered a fork (sidecar
`WeakKeyDictionary[Figure, RenderContext]` vs `fig._nmd_render_context`) and
recommended the sidecar. **Reject both**, for the same reasons, now applied to
`Style`:

- The codebase already passes the containing-width explicitly
  (`native_width_in=NATIVE_W` at every call site, e.g. SF33 `figure_sf33.py` main,
  SF42 `figure_s_model_comparison.py:266`). Attaching state to the figure
  *reverses* a working pattern.
- Figure-attached state keyed on a reused `fig` (composite panels are often built
  by mutating and re-validating a figure) is a stale-state bug waiting to happen:
  panel N's style lingers when panel N+1 is validated.
- Explicit params are trivially testable — a `Style` literal drives a rule with no
  matplotlib figure at all.

**Why `Style` is the right unifying abstraction — the identical-formula evidence.**
The four readability helpers in `ggplot_style.py` are two copies of the same two
functions, differing only in a width constant:

- `docx_body_fs` (:80) and `grant_body_fs` (:322) are `round(target_pt ×
  native_w / content_w)` with `content_w` = 6.5 vs 7.5.
- `docx_effective_pt` (:95) and `grant_effective_pt` (:337) are **the exact same
  formula** `native_pt × content_w / native_w`, again differing only in the width
  constant (6.5 vs 7.5).

Two media, one formula, one constant that varies: that is precisely a value that
belongs as a field of a `Style`, and the duplicated helpers collapse to
`style.body_fs(native_w)` and `style.effective_pt(native_pt)`. `PAPER` and `GRANT`
are the two style overlays that already exist as prose —
`~/.claude/memory/figures_style_publications.md` and `figures_style_grants.md` —
made executable.

**Design — a `Style` bundles a render envelope and a visual identity.**

```python
# figures/lib/style.py
@dataclass(frozen=True)
class RenderEnvelope:
    """Differs paper vs grant — the media-scaling spec. Formerly RenderContext."""
    content_w_in: float      # width the artifact ships at (docx 6.5 / NIH page 7.5)
    target_pt: float         # body target at content scale
    floor_pt: float          # readability floor at content scale
    header_pt: float         # header target at content scale

@dataclass(frozen=True)
class VisualIdentity:
    """Usually SHARED across a project's paper + grant."""
    palette: dict[str, str]  # semantic color map (NMD_COLOR, CONTROL_COLOR, CT_COLORS, …)
    body_font: str
    header_font: str
    spines: str              # spine/axis treatment (e.g. "ggplot_grey")
    panel_letter: PanelLetterSpec

@dataclass(frozen=True)
class Style:
    envelope: RenderEnvelope
    identity: VisualIdentity
    name: str

    def body_fs(self, native_w_in: float) -> int:
        return round(self.envelope.target_pt * native_w_in / self.envelope.content_w_in)

    def effective_pt(self, native_pt: float, native_w_in: float) -> float:
        return native_pt * self.envelope.content_w_in / native_w_in   # one formula, both media

    def with_(self, *, envelope=None, identity=None, name=None) -> "Style":
        return replace(self, envelope=envelope or self.envelope,
                       identity=identity or self.identity, name=name or self.name)
```

**Styles compose so the visual identity is not duplicated (drift, one level down).**
The envelope differs by medium but the visual identity is shared. Hardcoding
palette/fonts into both `PAPER` and `GRANT` would re-introduce exactly the kind of
copy-drift §4 is trying to eliminate — a palette edit that lands in one and not the
other. So define one project base style and override only the envelope per medium:

```python
PAPER_ENVELOPE = RenderEnvelope(content_w_in=PAPER_CONTENT_W_IN, target_pt=PAPER_TARGET_PT,   # :70-77
                                floor_pt=PAPER_FLOOR_PT, header_pt=HEADER_FS_AT_PAPER_SCALE)
GRANT_ENVELOPE = RenderEnvelope(content_w_in=GRANT_PAGE_W_IN, target_pt=GRANT_TARGET_PT,       # :316-319
                                floor_pt=GRANT_FLOOR_PT, header_pt=GRANT_HEADER_PT)

NMD_BASE = Style(envelope=PAPER_ENVELOPE,                       # envelope is a placeholder here
                 identity=VisualIdentity(palette={...NMD_COLOR, CONTROL_COLOR, CT_COLORS...},
                                         body_font="Arial", header_font="Arial",
                                         spines="ggplot_grey", panel_letter=PANEL_LETTER_SPEC),
                 name="nmd_base")
PAPER = NMD_BASE.with_(envelope=PAPER_ENVELOPE, name="PAPER")   # figures_style_publications.md, executable
GRANT = NMD_BASE.with_(envelope=GRANT_ENVELOPE, name="GRANT")   # figures_style_grants.md, executable
```

The envelope constants are the existing module constants verbatim:
`PAPER_CONTENT_W_IN=6.5`, `PAPER_TARGET_PT=10`, `PAPER_FLOOR_PT=7`
(ggplot_style:70–77); `GRANT_PAGE_W_IN=7.5`, `GRANT_TARGET_PT=11`,
`GRANT_FLOOR_PT=8`, `GRANT_HEADER_PT=13` (ggplot_style:316–319).

`render_and_validate` gains `style: Style = PAPER` and keeps `native_width_in`
accepted as a **deprecated fallback** during migration (when a caller passes only
`native_width_in`, reconstruct `PAPER` and read the width from the argument).
`assert_docx_readable`'s current `native_width_in`/`content_w_in`/`floor_pt` kwargs
(validate:1117–1132) become `style: Style | None = None` + `native_width_in`;
when `style is None` it defaults to `PAPER`, so no call site breaks. Composite
panels use `composite_effective_pt(panel_native_pt, panel_native_w, composite_w,
style)` (§3.3). Rules read `style.envelope.floor_pt` / `style.effective_pt(...)` /
`style.identity.palette`, never fig-attached state.

### 2.3 Verdict-based migration oracle (ACCEPTED — rewrites v1 Phase 0)

v1 Phase 0 (v1:55) captured "SHA-256 of each PDF and PNG, plus captured stdout."
The PDF SHA is unusable: matplotlib writes `CreationDate` and `ModDate` into the
PDF `/Info` dict and randomizes font-subset tags, so the hash changes on every
render regardless of pixels. Replace with:

- **Primary oracle — structured verdicts.** A harness renders all 17 gated SFs and
  records, per figure, a JSON record: for each check, `{pass|fail, [offender
  strings]}`. Two runs are equal iff every check on every figure has the same
  pass/fail and the same offender list. This is the invariant every *refactor*
  phase must hold, and the invariant every *fix* phase must hold **except for an
  explicitly enumerated set of expected flips** (Phase 1 is designed to flip
  specific figures pass→fail).
- **Secondary oracle — PNG pixels.** SHA-256 of the PNG (Agg output is
  deterministic; PNGs carry no date stamp unless one is passed). Any PNG hash
  change on a refactor phase is a real regression to investigate.
- **PDF comparison** only when needed, via `savefig(..., metadata={"CreationDate":
  None, "ModDate": None})` to strip the volatile keys, then compare; or rasterize
  the PDF and diff pixels. PDF hashes are never the gate.

Harness lives at `figures/lib/_verdict_baseline.py` (one-off, not part of the
gate); output at `figures/lib/_baseline/<sha>.json`.

### 2.4 Coverage-gap detection (ACCEPTED — new, not in v1)

`UnknownElement` (v1 §2.A) fires when the enumerator meets an artist class it has
no adapter for. It is silent when a *known* element pair or group is present that
**no rule matches** — e.g. a `PanelLetter` co-located with a `TickLabel` on the
same axis with no rule relating them (exactly the SF34 "C" overlapping ytick "60"
class the multipanel validator special-cases at validate:945–952). That hole is
invisible today.

Add `assert_registry_covers(fig, registry)`:

- Enumerate elements; compute (a) the set of element **types** present, (b) the set
  of unordered element-**type-pairs** that actually co-occur on the same
  owner-axis or within `fig`, (c) the set of group **roles** the registry's
  `GroupRule` selectors emit.
- For each present type with no `UnaryRule`, each present co-occurring type-pair
  with no `PairwiseRule`, and each nonempty group with no `GroupRule`, emit an
  `INFO`/`WARNING` "coverage gap: (PanelLetter, TickLabel) present on ax3, no
  rule applies."
- Run it in a repo-level report (not necessarily the per-figure gate) so the
  registry's blind spots are enumerable at any time.

This converts "we forgot to write a rule" from an invisible miss into a printed
line.

---

## 3. Target architecture

Four layers, bottom-up.

### 3.1 Typed enumerator (`figures/lib/elements.py`)

`enumerate_elements(fig) -> list[FigureElement]`, one pass over the figure, one
`fig.canvas.draw()`. Element classes map 1:1 to artists actually present across
SF25–SF42:

| Class | Backing artist | Today's ad-hoc source |
|---|---|---|
| `TickLabel` | x/y `get_xticklabels()` | multipanel Phase A (validate:864), style_symmetric (validate:1064) |
| `AxisLabel` | `ax.xaxis.label`/`yaxis.label` | validate:872, 1070 |
| `Title` | `ax.title` | validate:872 |
| `PanelLetter` | `ax.text` with `transform=ax.transAxes` from `panel_letter()` | is_axes_frac flag (validate:113) |
| `FacetHeader` | `ax.text` with `bbox=STRIP_BG` from `facet_header()` | validate:113 |
| `AxesText` | other `ax.texts` (sig stars, curve labels, medians) | validate:104 |
| `FigureText` | `fig.texts` / suptitle | validate:890 |
| `Legend` | **`ax.get_legend()` AND `fig.legends`** | validate:1240, 1597 (fig.legends missing today) |
| `DataMask` | rasterized data pixels for an axis | `_compute_data_mask` (validate:1214) |
| `Shape` (`kind ∈ rect/polygon/circle`) | patches | validate:196 |
| `Segment` | arrows/Line2D | validate:151 |
| `UnknownElement` | anything unmatched | — |

Every element exposes `owner_ax`, `bbox_px` (display coords via a shared
`BboxPx.from_display`), `bbox_data`, and any role tags (`is_axes_frac`,
`kind`). Element `.integrity()` hosts unary self-checks. The single most important
correctness win: `Legend` enumerates **both** `ax.get_legend()` and each
`fig.legends` entry — closing SF42 (§5) structurally rather than by a one-off loop.

### 3.2 N-ary-capable rule layer (`figures/lib/rules/`)

The three rule kinds from §2.1, in `rules/model.py`; concrete rules in
`rules/layout.py`, `rules/readability.py`, `rules/style.py`. A `Registry`
aggregates `UnaryRule`/`PairwiseRule`/`GroupRule` instances and answers coverage
queries (§2.4).

### 3.3 Explicit `Style` (`figures/lib/style.py`)

Per §2.2. A `Style` (render envelope + visual identity) is passed into every
rule's `fn`. Readability rules read `style.envelope.floor_pt`,
`style.effective_pt(native_pt, native_w_in)`; style/identity rules read
`style.identity.palette` / `style.identity.body_font`. Composite panels use
`composite_effective_pt(panel_native_pt, panel_native_w, composite_native_w,
style)` — the panel's own `native_width_in` lies, so the effective size is
computed against the composite's shipping width using `style`'s envelope:

```python
def composite_effective_pt(panel_native_pt, panel_native_w, composite_native_w, style):
    # panel drawn at panel_native_w, tiled into a composite that ships at
    # style.envelope.content_w_in scaled from composite_native_w.
    return panel_native_pt * (panel_native_w / composite_native_w) \
           * (style.envelope.content_w_in / composite_native_w)
```

`PAPER` and `GRANT` are the two shipped `Style` instances (§2.2); the type is open
for future media.

### 3.4 Runner + gate

```python
def run_rules(fig, style, registry) -> list[Issue]:
    els = enumerate_elements(fig)
    issues = []
    for e in els:
        for r in registry.unary_for(type(e)):
            if (i := r.fn(e, style)): issues.append(i)
    by_pair = _cooccurring_pairs(els)          # same owner_ax or same fig
    for (a, b) in by_pair:
        for r in registry.pairwise_for(type(a), type(b)):
            if (i := r.fn(a, b, style)): issues.append(i)
    for r in registry.group_rules:
        for group in r.selector(els):
            issues += r.fn(group, style)
    issues += assert_registry_covers(els, registry)   # coverage, §2.4
    return issues
```

`render_and_validate` (ggplot_style:100) gains `style: Style = PAPER` (with
`native_width_in` kept as a deprecated fallback, §2.2); internally it calls
`run_rules(fig, style, registry)`. Each public `assert_*` becomes a thin shim that
runs the single corresponding rule and raises on any ERROR — so every SF import
(`from ggplot_style import render_and_validate`, SF42:41) and every direct
`assert_*` call keeps working.

### 3.5 Worked conversions

#### (a) N-ARY: `assert_style_symmetric` → `GroupRule`

Today (validate:1040–1100): iterate `fig.axes`, build `per_axis[i][role]`, then
per role compare `max-min` against tolerance. As a group rule the selector yields
one group per role (each group = the `AxisFont` readings for that role across all
sibling axes), and `fn` judges the spread of the whole group:

```python
# rules/style.py
def _by_role(els):
    fonts = [e for e in els if isinstance(e, AxisFont)]   # AxisFont: (ax_idx, role, size)
    groups = defaultdict(list)
    for f in fonts:
        groups[f.role].append(f)                          # role ∈ xtick/ytick/xlabel/ylabel
    return [g for g in groups.values() if len(g) >= 2]

def _style_symmetric(group, style, *, tol_pt=0.5):
    sizes = sorted({round(f.size, 1) for f in group})
    if sizes[-1] - sizes[0] <= tol_pt:
        return []
    role = group[0].role
    parts = " vs ".join(f"{s}pt on axes {[f.ax_idx for f in group if round(f.size,1)==s]}"
                        for s in sizes)
    return [Issue("style_symmetric", "ERROR",
                  f"{role}: {parts}",
                  "Every panel must use the same BODY_FS from ggplot_style; "
                  "remove local literals (fontsize=11, BODY_FS-2, …)")]

STYLE_SYMMETRIC = GroupRule(selector=_by_role, fn=_style_symmetric, name="style_symmetric")
```

Note this is the natural home: the group *is* the sibling set, and the report
lists which axes hold which size — impossible to express as a pairwise rule
without O(n²) redundant emissions. The alignment check (validate:592) converts the
same way (selector = greedy y-band partition, `fn` = deviation-from-band-mean over
groups of ≥3).

#### (b) PAIRWISE: `assert_legend_clear` → `PairwiseRule(Legend, DataMask)`

Today (validate:1570–1677): `leg = ax.get_legend()` (validate:1597 — **misses
`fig.legends`**), compute mask, dilate, convert legend bbox to mask coords, count
intersection, build a slot suggestion. As a pairwise rule the pixel arithmetic
lives on the elements and the rule is the geometry:

```python
# rules/layout.py
def _legend_clearance(legend: Legend, data: DataMask, style, *,
                      clearance_px=4, overlap_tol_px=0):
    if legend.owner_ax is not data.ax:      # only clear against its own axis' data
        return None
    n = data.dilated(clearance_px).intersect_px(legend.bbox_px)
    if n <= overlap_tol_px:
        return None
    return Issue("legend_clearance", "ERROR",
                 f"Legend overlaps data (or within {clearance_px}px): {n}px",
                 data.suggest_open_slot(legend.bbox_px))   # was 40 lines inline

LEGEND_CLEARANCE = PairwiseRule(Legend, DataMask, _legend_clearance, symmetric=False)
```

Because `Legend` now enumerates `fig.legends` too (§3.1), and `_compute_data_mask`
hides all `fig.legends` before rasterizing (Phase-1 fix, §5), SF42's figure-level
legend is checked against every axis' data mask automatically — the defect closes
in the enumerator, not in the rule.

---

## 4. Packaging / single source of truth

### 4.1 The drift is already real and already broken

There are two copies of the toolchain and they have diverged badly:

- `~/.claude/utils/validate_figure_layout.py` — **34,854 bytes, dated 2026-06-13.**
  Its last function is `validate_figure_layout` (line 401). It is **missing**
  `validate_multipanel_layout`, `assert_style_symmetric`, `assert_docx_readable`,
  `_compute_data_mask`, `assert_data_free`, `find_open_regions`,
  `assert_tick_labels_disjoint`, and `assert_legend_clear` — i.e. all eight of the
  functions the current gate depends on.
- `figures/lib/validate_figure_layout.py` — **67,130 bytes, dated 2026-07-10.**
  Canonical; nearly 2× the size.

Worse, the split is *complementary and incoherent*: `figure_lint.py` and
`figure_template.py` exist only under `~/.claude/utils/`; `ggplot_style.py` (which
defines `render_and_validate`, `docx_body_fs` — the very symbols `figure_lint.py`
greps for) exists only under `figures/lib/`. Neither copy is a working whole. The
per-repo-copy model has already guaranteed the drift the critique warned about.
CLAUDE.md already tells authors "use `figures/lib/` (not `~/.claude/`)", but grant
figures live in `~/claude_projects/grants/nmd_2026/` and import `ggplot_style`'s
grant helpers — so *something* has to resolve the library from outside this repo.

### 4.2 Options

| Option | SSOT? | Reproducible for Yul (fresh clone)? | Resolves cross-repo (grants)? | Drift risk |
|---|---|---|---|---|
| **A. Vendored per-repo copy** (status quo) | No | Yes | No (each repo re-copies) | **Proven high** |
| **B. Symlink** `figures/lib` → canonical checkout | Yes | **No** (symlink dangles on clone) | Fragile | Low but breaks repro |
| **C. Installable package** `nmd-figtools`, version-pinned | **Yes** | Yes (pinned in env) | **Yes** | Low |

### 4.3 Recommendation — Option C

Extract `figures/lib/` (validator + `ggplot_style` + `figure_geometry` +
primitives + `figure_lint`) into a standalone installable package
**`nmd-figtools`**, in its own repo, **dual-pushed to public GitHub + Channing
GitLab** — mirroring the pattern the project already uses for `Isopair` and
`NMD_orf_model_v5_4ct` (per this repo's CLAUDE.md "Linked repos"). Consuming
repos declare a pinned version:

- The NMD paper repo installs `nmd-figtools==X.Y` into its conda env; `figures/lib`
  becomes a thin re-export shim (`from nmd_figtools import *`) OR is removed and
  imports point at the package. Yul gets a reproducible figure toolchain by
  `pip install -r requirements.txt`, version pinned.
- The grants tree installs the same pinned package; grant helpers resolve from one
  place.
- **The package must export the `Style` instances themselves** — `PAPER` and
  `GRANT` (and the `NMD_BASE` they compose from, §2.2) — so both the paper repo
  (`nmd/figures`) and the grant repo (`grants/nmd_2026`) do `from nmd_figtools
  import PAPER, GRANT` and validate against the *same* palette, typography, and
  envelope objects. This is the composition point where drift would otherwise
  reappear: if each repo defined its own `PAPER`, the palette would fork exactly as
  the validator copies did. One exported `Style` per medium is the SSOT for visual
  identity, not just for validator code.
- **Delete `~/.claude/utils/validate_figure_layout.py`** (it is a stale partial and
  actively dangerous — importing it silently loses eight validators). If ambient
  figure work outside any repo still needs the tools, make `~/.claude/utils/` a
  re-export of the installed package, never a fork.

Add a **drift guard** for the transition period: a `nmd_figtools.__version__` and a
`figures/lib/_sync_check.py` that asserts the installed package version matches the
version pinned in the repo; wire it into the same lint pass as `figure_lint.py` so
CI fails on drift. This is the "vendored-with-sync-check" fallback *only* if a
fully-installable flow is delayed; the target is a clean `pip install`.

Tradeoff acknowledged: a package adds a release step (bump version, dual-push, pin
in consumers). That cost is trivial against the demonstrated failure mode of eight
silently-missing validators, and it is the *only* option that is simultaneously
SSOT, clone-reproducible, and cross-repo.

---

## 5. Phased migration

**Phase 1 is architecture-agnostic** — every edit is a surgical change to an
existing function that is correct under any future design, and each closes a named
open defect. Phases 2–6 are substrate and are verdict-preserving by construction
(§2.3), except where a phase deliberately adds a new check.

### Phase 1 — In-place defect fixes (architecture-agnostic)

Correct under v1's design, v2's design, or none. Do these first; they retire the
open figure defects and give the substrate work a richer verdict baseline to
preserve.

1. **Extend `assert_text_within_canvas` walk** (ggplot_style:472–473). Today it
   walks only `ax.texts + [ax.xaxis.label, ax.yaxis.label, ax.title]`. Add: tick
   labels (`get_xticklabels()`/`get_yticklabels()`), `ax.get_legend()`,
   `fig.legends`, and `fig.texts`. → **Closes SF31 legend truncation** (legend
   text clipped at the canvas edge is currently invisible to the canvas guard).
2. **Loop `fig.legends` in the legend + data-mask checks.**
   - `assert_legend_clear` (validate:1597): iterate `ax.get_legend()` *and* every
     `fig.legends` entry. → **Closes SF42 figure-level-legend miss**
     (`figure_s_model_comparison.py:242` uses `fig.legend(...)`, never checked).
   - `_compute_data_mask` hide-list (validate:1240): also hide all `fig.legends`
     before rasterizing. → **Closes the legend-as-data mask** (a figure-level
     legend's pixels are currently counted as "data," corrupting both
     `assert_data_free` and `assert_legend_clear` for every axis beneath it).
3. **Broaden `assert_tick_labels_disjoint` past adjacent-only** (validate:1533,
   `for i in range(len(labels) - 1)`). A wide rotated label can overlap the label
   two positions away; adjacent-only misses it. Switch to all same-axis pairs
   (n is small per axis) or a sweep-line. → **Closes the SF29 / SF33 / SF34
   tick-crowding gaps.**
4. **Re-tune the `assert_data_free` skip heuristic** (validate:1363–1366). Today
   `sub.size/axes_area < 0.015` **or** `len(s.strip()) <= 8` → skip. The short-
   string clause is too permissive (it skips short annotations that genuinely sit
   on data) and interacts with the legend-as-data bug. Make the short-string skip
   conditional on *also* being small-bbox and *not* part of a legend; keep the
   small-bbox skip for true value markers. → **Closes the B3 legend-overlap
   residuals (SF36 / SF38 / SF40 / SF41).**

**Oracle.** Run the §2.3 verdict harness before and after. Phase 1 is *expected* to
flip a specific, enumerated set of figures pass→fail (SF31, SF42, SF29/33/34,
SF36/38/40/41). Write that expected-flip list first; any figure changing verdict
that is *not* on the list is a regression. Then fix each newly-failing figure on
its own commit (this is the point of the phase). **Rollback:** each fix is one
function; revert individually.

### Phase 2 — Package extraction / SSOT (§4)

Extract `figures/lib` into `nmd-figtools`, pin it in the paper + grants repos, make
`figures/lib` a re-export shim, delete/redirect the stale `~/.claude/utils` copy,
add the sync-check guard. **Oracle:** verdicts + PNG hashes identical across all 17
SFs before/after the import switch. **Rollback:** repin to the vendored copy.

### Phase 3 — Typed enumerator + coverage report (additive)

New `elements.py` (§3.1) and `assert_registry_covers` (§2.4). Nothing wired into
the gate yet. **Oracle:** enumerate all 17 SFs → zero `UnknownElement`; produce and
review the coverage-gap report. **Rollback:** delete the files.

### Phase 4 — Explicit `Style` (additive, backward-compatible)

`style.py` (§2.2, §3.3): `RenderEnvelope` + `VisualIdentity` + `Style`, the
`NMD_BASE`-composed `PAPER`/`GRANT` instances, and the collapse of
`docx_body_fs`/`grant_body_fs` → `style.body_fs(...)` and
`docx_effective_pt`/`grant_effective_pt` → `style.effective_pt(...)` (same formula,
§2.2). `render_and_validate` gains `style: Style = PAPER`; `assert_docx_readable`
gains `style: Style | None = None` defaulting to `PAPER`; `native_width_in` kept as
a deprecated fallback so no call site changes. Add `composite_effective_pt(...,
style)`. **Oracle:** verdicts identical (the default `PAPER` path reconstructs
today's behavior). **Rollback:** drop the `style` param; keep the old helpers.

### Phase 5 — Rule model + registry + runner; convert `assert_*` to shims

`rules/` package (§3.2, §3.4). Convert each validator to the corresponding
Unary/Pairwise/Group rule (worked examples §3.5); public `assert_*` names become
thin shims; wire `assert_registry_covers` into the repo-level report.
**Oracle:** verdicts identical **except** the intentional coverage-gap
INFO/WARNING lines, which are enumerated up front. This is the phase most at risk
of subtle verdict drift — run the harness per converted check, not just at the
end. **Rollback:** the shims are ~10 lines each; revert per check.

### Phase 6 — New checks the substrate unlocks (opt-in)

- `assert_composite_font_effective` using `composite_effective_pt` — catches the
  SF33-family bug where a panel embedded in a composite is measured against its own
  `native_width_in` instead of the composite's shipping width.
- `assert_palette_consistent(figs)` — cross-figure `GroupRule` over `FigureText`/
  `AxesText` colors, run from a repo-level runner.

Opt-in; not in the per-figure gate until each is validated against all 17 SFs.
**Oracle:** each new check re-runs the baseline; surprises triaged before mainline.

**Sequencing rationale.** Phase 1 pays down the open defects immediately and
enriches the verdict baseline (more real failures captured) *before* the refactor,
so Phases 2–5 have a stricter invariant to preserve. Every phase is independently
revertable; no phase requires the next. `render_and_validate`'s signature
(ggplot_style:100) is unchanged throughout.

---

## 6. Risks & open questions

- **Phase 5 verdict drift.** Moving pixel arithmetic from inlined validators onto
  element methods can change rounding at bbox boundaries (e.g. the
  `floor`/`ceil` conversions at validate:1617–1620). Mitigate: port the exact
  integer math into `BboxPx`/`DataMask` unchanged; run the verdict harness per
  converted check.
- **`fig.legends` bbox in `_compute_data_mask`.** A figure-level legend spans
  multiple axes; hiding it is correct, but its clearance must then be checked
  against *each* underlying axis' mask (the `owner_ax` guard in `_legend_clearance`
  must treat a `fig`-owned legend as owned-by-all). Decide the ownership semantics
  in Phase 1 (the fix) and carry them into Phase 5 (the rule).
- **Package boundary for `ggplot_style`.** `render_and_validate` lives in
  `ggplot_style` and imports the validators lazily to avoid a cycle
  (ggplot_style:157–166). In `nmd-figtools` decide whether `render_and_validate`
  stays in the style module or moves to a `gate` module; either works, but pick
  before Phase 2 to avoid a second import rearrangement.
- **Coverage-report noise.** `assert_registry_covers` will initially flag many
  benign type-pairs (every `PanelLetter`×`TickLabel` that doesn't actually
  collide). Decide whether coverage is INFO-only (report) or WARNING (nag), and
  whether to seed an explicit "intentionally uncovered" whitelist — mirror the
  same-owner whitelist already encoded at validate:945–952.
- **PNG determinism assumption.** The secondary oracle assumes Agg PNG output is
  bit-stable across the same matplotlib/freetype versions. If Yul's environment
  differs, PNG hashes will diff for environment reasons, not code reasons — pin
  matplotlib/freetype in `nmd-figtools`' requirements, or fall back to
  tolerance-based pixel diff.
- **Grant vs paper floor divergence.** `PAPER_FLOOR_PT` was lowered 9→7 on
  2026-07-10 (ggplot_style:77) for tilted violin ticks; `GRANT_FLOOR_PT` is 8. The
  `PAPER`/`GRANT` `Style` envelopes cleanly carry both, but confirm no SF silently
  relies on the old 9pt floor before Phase 4 folds the constant into the envelope.
- **Style-composition boundary — what is shared vs per-medium.** `NMD_BASE` assumes
  palette + typography + spines are identical paper-vs-grant. Confirm this against
  the two prose overlays (`figures_style_publications.md` vs
  `figures_style_grants.md`) — if any identity attribute legitimately differs by
  medium (e.g. grant panel-letter sizing), it must move from `VisualIdentity` into
  the per-medium override, not be forced shared. Decide the split before Phase 4.
