# Scientific-Figure Validator: Architecture Plan (v1 — Plan Agent 1)

## 1. Executive Summary

**Proposal.** Extract a typed element enumerator (`figures/lib/elements.py`) that walks the matplotlib figure once and returns a heterogeneous list of `FigureElement` instances (one subclass per real artist class in your 15 SFs). Each element carries its own integrity checks; a small pairwise rule registry expresses relationship checks (overlap / clearance / alignment) as `(element_type_a, element_type_b) → callable`. Readability collapses to a `RenderContext` value object attached to the figure that encodes the containing space (paper 6.5", grant 7.5", composite panel-of-a-composite) — every readability floor derives from that context rather than a per-figure `native_width_in` argument. `render_and_validate` stays a one-liner; new element classes are added in one file; new relationship checks are one decorator away and automatically cover every element pair the registry says they apply to.

**Why.** Every one of the five recent blind spots — the SF31 legend truncation, the SF42 opt-in miss, the legend-as-data mask, the tick-crowding gap, the density-tail overlap — shares one shape: a rule that lives inside one validator's private artist walk. Fix the walk (one enumerator, typed elements, a place a new element class has to be declared) and every existing rule inherits coverage. Fix the readability model (containing-space value object) and multipanel composites stop needing a per-panel `native_width_in` that lies. Fix the tiering (single / cross-figure / composite) by making tier-2 and tier-3 explicit callable groups instead of hoping each figure's `main()` remembers to call them. None of this requires touching the `render_and_validate` signature.

## 2. Three Candidate Architectures

### A. Element Registry + Typed Enumerator + Pairwise Rule Table

**Model.** One module (`elements.py`) with `class FigureElement` and subclasses that map 1:1 to artists actually present in the SFs today: `TickLabel`, `AxisLabel`, `Title`, `PanelLetter` (via `panel_letter()`), `AxesText` (`ax.texts` that isn't a panel letter — SF29 significance stars, SF33 curve labels, violin median labels), `FacetHeader` (`ax.text` with `bbox=STRIP_BG`, produced by `facet_header()` in SF33-family), `Legend`, `FigureText` (`fig.texts` / suptitle), `DataMask` (the rasterized data-drawn mask), and the four schematic-shape classes already there (`Rect`, `Segment`, `Polygon`, `Circle`). Each element has an `.integrity()` method returning `list[Issue]`. A separate `RULE_REGISTRY: dict[frozenset[type], list[Rule]]` holds relationship rules; `run_checks(fig)` walks all pairs and dispatches. Unknown artist types show up as `UnknownElement` and produce a WARNING — a new artist class cannot silently pass.

**Addresses direction (1).** Two-axes checks are structurally enforced: `.integrity()` is per-element; the rule table is pairwise. **Direction (2):** `RenderContext` (a dataclass with `content_w_in`, `role`, `floor_pt`, `target_pt`) lives on `fig` (e.g. `fig._render_context`) and is queried by every readability rule. Composite panels inherit their composite's context automatically. **Direction (3):** three explicit runner groups — `run_single(fig)`, `run_composite(fig)`, `run_cross_figure(figs)`.

**Tradeoffs.** Adoption cost: moderate; a shim in the current `assert_*` functions preserves signatures during migration. Blind-spot resistance: **high** — new element classes must register, and `UnknownElement` catches accidental omissions. Ergonomics: `render_and_validate` unchanged. Testability: any rule is a pure function of typed elements; a fake `Legend(bbox=…)` fixture drives it without matplotlib.

### B. Rule-DSL / Declarative Constraints

**Model.** Rules declared as strings/decorators: `@rule("no_overlap(Legend, DataMask, clearance=4px)")`. A parser resolves them to element-set matches. Analogous to how CSS selectors match DOM elements. `assert_legend_clear` becomes one line of DSL, and `render_and_validate` runs a rule file.

**Addresses (1):** cleanly, at the cost of parsing. **(2):** context expressions like `at(context.docx)` need a mini-eval. **(3):** tiers become rule scopes: `@scope("single")`, `@scope("composite")`.

**Tradeoffs.** Elegant for adding new checks. Debuggability suffers — stack traces point into the DSL engine, not user code. Adoption cost high (rewrite all seven checks + build/maintain the DSL). Ergonomics for figure authors unchanged. Testability suffers: you now test the DSL runtime as well as the rules. **Rejected** as over-engineered for a codebase whose largest asset is Python.

### C. Incremental Refactor with Centralized Enumeration Only

**Model.** Keep every existing `assert_*` function. Add one helper `enumerate_elements(fig) → list[TypedElement]`. Rewrite each validator internally to consume the enumerated list instead of walking artists itself. No registry, no rule table, no context object.

**Addresses (1):** integrity checks moved into enumeration; relationship checks remain per-validator. **(2):** no built-in mechanism; each validator still takes `native_width_in`. **(3):** no formal tiering; you'd still call three groups by convention.

**Tradeoffs.** Adoption cost: **lowest.** Blind-spot resistance: **medium** — enumeration catches "which artists exist" blind spots (fixing the SF31 legend miss centrally); still lets a new pairwise blind spot slip if a rule forgets to iterate the new element. Ergonomics: unchanged. Testability: unchanged. **Half-measure:** it solves the enumeration blind-spot class but leaves the pairwise blind-spot class open.

## 3. Recommendation

**Architecture A, with two concessions borrowed from C: (i) migrate incrementally rather than big-bang; (ii) keep `assert_*` public names as thin wrappers around the registry so `render_and_validate` and every SF's own imports don't change.**

Rationale:
- Every real blind spot you listed is a walk/enumeration bug. A typed enumerator with an `UnknownElement` fallback is the smallest thing that structurally prevents them.
- The `RenderContext` isn't optional for direction (2). Threading `native_width_in` through composite panels breaks — a panel inside an 11" composite that itself will be scaled to 6.5" docx needs *composite* context, not its own. This is only clean with a context object attached to the figure.
- The rule registry is small (~15 rules today, ~40 at maturity) and dispatch is O(n²) over element pairs, which is fine — the SF33 figure has ~30 elements, so ~450 pair checks, all cheap once bboxes are cached.
- DSL is deferred until the rule table exceeds ~50 entries or you start writing python that reads like YAML anyway. That day may never come.

**Fork you should decide before Phase 3.** How to attach `RenderContext` to `fig`. Two options:
- **(a) Sidecar dict:** `_RENDER_CTX: WeakKeyDictionary[Figure, RenderContext]` maintained by a `set_render_context(fig, ctx)` helper. Zero attribute pollution on matplotlib objects.
- **(b) Attribute:** `fig._nmd_render_context = ...`. Simpler, but attaches to matplotlib's Figure and can confuse other tooling that introspects `fig`.

Recommend (a); flag both in the plan and pick before Phase 3 starts.

## 4. Phased Migration Plan

Baseline invariant across all phases: the 15 currently-migrated SFs must produce byte-identical PDFs and PNGs (or the diff must be explicitly reviewed and accepted). Baseline captured in Phase 0.

**Phase 0 — Baseline (0.5 day).** Add `figures/lib/_snapshot_baseline.py` (a one-off script, not committed to the gate) that renders all 15 SFs and records: SHA-256 of each PDF and PNG, plus captured stdout of `render_and_validate`. This is the rollback oracle. **Rollback point:** any subsequent phase that fails snapshot comparison reverts to prior commit. **Detection:** snapshot diff.

**Phase 1 — Element enumerator (2 days).** Additive only. New file `figures/lib/elements.py`:
- `FigureElement` ABC with `.bbox_px`, `.bbox_data`, `.owner_ax`, `.integrity()`.
- Concrete subclasses for the 10 element classes above (all backed by artists in the 15 SFs).
- `enumerate_elements(fig) → list[FigureElement]` — one pass; `UnknownElement` for anything not matching a registered adapter.
- Unit-adjacent smoke test: run enumerator on each of the 15 SFs, assert no `UnknownElement`. Runs in CI/local before merging.

Nothing wired into `render_and_validate` yet. **Rollback:** delete the file. **Detection:** the smoke test.

**Phase 2 — RenderContext (1 day).** Additive to `ggplot_style.py`:
- `@dataclass class RenderContext: content_w_in, role, target_pt, floor_pt`.
- Predefined: `PAPER_DOCX_CONTEXT`, `GRANT_PAGE_CONTEXT`, and a factory `composite_panel_context(composite_ctx, panel_share)`.
- `set_render_context(fig, ctx)` + `get_render_context(fig)` using WeakKeyDictionary (option (a) above, pending fork decision).
- Old `PAPER_*`, `GRANT_*` constants remain and become fields of the two contexts. `docx_body_fs`, `docx_effective_pt`, `grant_body_fs` unchanged.

No SF touches this yet. **Rollback:** delete the additions; constants are unchanged. **Detection:** none needed — no behaviour changed.

**Phase 3 — Convert `assert_*` to consume enumerator (2 days).** For each of `assert_text_within_canvas`, `assert_docx_readable`, `assert_data_free`, `assert_legend_clear`, `assert_tick_labels_disjoint`, `assert_style_symmetric`, plus the two `validate_*_layout` monoliths: rewrite the internals to walk `enumerate_elements(fig)` instead of hand-rolled artist iteration. Signatures unchanged. `assert_docx_readable` gains an internal fallback: if `get_render_context(fig)` is set, use it; else fall back to `native_width_in` — that's how you preserve the current SF wiring without changing any SF.

After the rewrite, re-render all 15 SFs and diff against the Phase 0 baseline. Any diff must be reviewed (many will be no-ops; some may show tighter clearance because the enumerator now finds an artist a previous walk missed — SF31 is the leading candidate). **Rollback:** revert `validate_figure_layout.py`. **Detection:** baseline diff.

**What breaks and how you spot it:** the SF31 legend-truncation bug will now fire, because `assert_text_within_canvas` will pick up `Legend` via the enumerator. That's a *feature* — but it means SF31 fails re-render until fixed. Any figure the enumerator finds an artist for that was previously invisible to a validator will similarly regress. Detected by baseline snapshot script. Fix each regression on its own commit before moving to Phase 4.

**Phase 4 — Rule registry + tiered runner (1 day).** Introduce `RULE_REGISTRY` and rewrite `render_and_validate` to call one of three runners: `run_single_figure(fig)`, `run_composite(fig)`, `run_cross_figure(figs)`. `render_and_validate` still calls `run_single_figure` + `run_composite`. Old kwargs (`run_data_free`, `run_legend_clear`, `run_tick_disjoint`) become rule enable/disable flags on the runner. **Rollback:** the runners are 20 lines each; revert `ggplot_style.render_and_validate`. **Detection:** baseline diff.

**Phase 5 — New checks that the tiering unlocks (1 day, opt-in).**
- `assert_composite_font_effective` — computes each panel's effective docx pt from the *composite's* RenderContext, not the panel's own `native_width_in`. Catches the SF33-family bug where a panel embedded in a composite is measured against the composite's shipping width.
- `assert_palette_consistent(figs)` — cross-figure tier-2 check that walks the FigureText/AxesText color of every SF in a group and flags off-palette colors. Called from a repo-level runner, not per-figure.
- `assert_panel_letter_style(fig)` — enforces the `PANEL_LETTER_FS` and `x, y` position convention across every panel in a composite.

These are opt-in: not called from `render_and_validate` until they've been validated against all 15 SFs. **Rollback:** stop calling them. **Detection:** each new check re-runs the 15-SF baseline; any surprises get triaged before the check goes into the mainline.

**Phase 6 — Cleanup (0.5 day).** Remove transitional shims (any `if get_render_context(fig) is None: fallback` branches whose fallback is dead). Delete `_snapshot_baseline.py`.

**Total: ~7 working days.** Each phase is independently revertable; no phase requires the next.

**Explicit non-breakage:** `render_and_validate`'s signature is unchanged through all six phases. Every SF's `main()` continues to work with zero edits. The only visible surface change to a figure author is that new checks (Phase 5) *may* start failing for real defects — which is what you want.

## 5. Worked Example: `assert_legend_clear` Before/After

The current implementation lives at `/Users/petecastaldi/claude_projects/nmd/figures/lib/validate_figure_layout.py` lines 1570–1677, roughly 110 lines. It hand-rolls: legend fetch + visibility check, data-mask compute, dilate, legend bbox in display coords, coord conversion to mask coords, area check, `find_open_regions` call, per-slot anchor-map lookup, message construction with the suggested `bbox_to_anchor`. Every unrelated concern is inlined.

**Before (essential shape, condensed):**

```python
def assert_legend_clear(fig, ax, *, clearance_px=4, overlap_tolerance_px=0):
    leg = ax.get_legend()
    if leg is None or not leg.get_visible():
        return
    mask, (mx0, my0, mx1, my1) = _compute_data_mask(fig, ax)
    if mask.size == 0:
        return
    dilated = _dilate_mask(mask, clearance_px)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lb = leg.get_window_extent(renderer=renderer)
    if lb.width <= 0 or lb.height <= 0:
        return
    fig_h = int(round(fig.bbox.height))
    lx0 = max(0, int(math.floor(lb.x0 - mx0)))
    lx1 = min(mask.shape[1], int(math.ceil(lb.x1 - mx0)))
    ly1 = max(0, int(math.floor(fig_h - lb.y1 - my0)))
    ly0 = min(mask.shape[0], int(math.ceil(fig_h - lb.y0 - my0)))
    if lx1 <= lx0 or ly0 <= ly1:
        return
    sub = dilated[ly1:ly0, lx0:lx1]
    n_overlap = int(sub.sum())
    if n_overlap <= overlap_tolerance_px:
        return
    # ... 40 lines of message construction with slot suggestion ...
    raise AssertionError(...)
```

**After.** The pixel arithmetic and slot-suggestion messaging move into shared infrastructure. The rule itself is the geometry check:

```python
# figures/lib/elements.py (excerpt)
@dataclass
class Legend(FigureElement):
    artist: matplotlib.legend.Legend
    owner_ax: Axes

    @property
    def bbox_px(self) -> BboxPx:
        return BboxPx.from_display(self.artist.get_window_extent(_renderer(self.artist.figure)))

    def integrity(self) -> list[Issue]:
        return _within_canvas(self)  # inherits the standard bbox-inside-canvas rule

@dataclass
class DataMask(FigureElement):
    ax: Axes
    _mask: np.ndarray
    _origin_px: tuple[int, int, int, int]

    def dilated(self, px: int) -> "DataMask":
        return replace(self, _mask=_dilate(self._mask, px))

    def intersect_px(self, bbox: BboxPx) -> int:
        return _count_intersection(self._mask, self._origin_px, bbox)

    def suggest_open_slot(self, for_bbox: BboxPx) -> str:
        return _format_slot_suggestion(find_open_regions_from_mask(self._mask), self.ax)
```

```python
# figures/lib/rules_layout.py
@rule(Legend, DataMask, tier="single", name="legend_clearance")
def legend_clearance(legend: Legend, data: DataMask, *,
                     clearance_px: int = 4,
                     overlap_tolerance_px: int = 0) -> Issue | None:
    if legend.owner_ax is not data.ax:
        return None
    n = data.dilated(clearance_px).intersect_px(legend.bbox_px)
    if n <= overlap_tolerance_px:
        return None
    return Issue(
        check="legend_clearance",
        severity="ERROR",
        message=f"Legend overlaps data (or comes within {clearance_px} px): {n} px",
        suggestion=data.suggest_open_slot(legend.bbox_px),
    )
```

```python
# figures/lib/validate_figure_layout.py — public shim stays put
def assert_legend_clear(fig, ax, *, clearance_px=4, overlap_tolerance_px=0):
    elements = enumerate_elements(fig)
    legend = _first(elements, Legend, owner_ax=ax)
    data = _first(elements, DataMask, ax=ax)
    if legend is None or data is None:
        return
    issue = legend_clearance(legend, data,
                             clearance_px=clearance_px,
                             overlap_tolerance_px=overlap_tolerance_px)
    if issue:
        raise AssertionError(issue.format())
```

**What improved concretely:**
1. The pixel-space conversion (`lb → lx0/lx1/ly0/ly1`) is written once, in `BboxPx.from_display`, used by every rule that needs a display-coord bbox.
2. The mask dilation is a method on `DataMask` — the tick-crowding check (which also wants a "keep-out zone around every X" concept) uses the exact same primitive.
3. The slot-suggestion formatting (currently 40 lines of `anchor_map` munging) becomes `data.suggest_open_slot(bbox)` — reusable by `assert_data_free` for its "move to open region" advice.
4. The rule itself is 8 lines and has type-checked inputs; a fake `Legend` and `DataMask` in tests drive it in isolation from matplotlib.
5. Adding "legend must clear other panels' data in a composite" is a new `@rule(Legend, DataMask)` variant scoped to `tier="composite"` — no new mask code, no new bbox code.
6. `assert_legend_clear`'s public signature is unchanged. `render_and_validate` doesn't move.

The refactored version is shorter, moves plumbing to reusable primitives, and the *rule content* — the two lines that say "dilate the mask, check the intersection count" — is now visible without scrolling past bbox conversion boilerplate. That's the credibility test: the check reads like the sentence you'd say out loud describing what it does.

## Critical Files for Implementation
- /Users/petecastaldi/claude_projects/nmd/figures/lib/validate_figure_layout.py
- /Users/petecastaldi/claude_projects/nmd/figures/lib/ggplot_style.py
- /Users/petecastaldi/claude_projects/nmd/figures/lib/elements.py (to be created)
- /Users/petecastaldi/claude_projects/nmd/figures/lib/rules_layout.py (to be created)
- /Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures/SF33_TD2Bias_broad/figure_sf33.py
