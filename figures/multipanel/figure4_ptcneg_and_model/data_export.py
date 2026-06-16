"""
Figure 4 — DL-model side data export (panels D–H).

Reads the per-window-size metrics JSONs + best-model DeepSHAP TSVs from the
NMD_orf_model_v5_4ct results directory, writes panel-ready TSVs to ./data/.

Source of truth for the underlying model artifacts:
    ~/claude_projects/nmd/results/isoform_transitions/NMD_orf_model_v5/results/

(That directory is the local rsync mirror of peter4244/NMD_orf_model_v5_4ct
results/. See CLAUDE.md → "Linked repos".)

Panels:
    D — ROC curve of the best model (mirrors orf_model_report fig1d_roc_curve)
    E — Branch-level Shapley (Structural / Stop / ATG) — mirrors fig2a_branch_importance
    F — ORF-structural feature Shapley — mirrors fig2b_structural_importance
    G — Nucleotide Shapley around START window — mirrors fig2c_atg_profile
    H — Nucleotide Shapley around STOP window — mirrors fig2d_stop_profile

Run from the figure folder:
    python data_export.py
"""

import json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "data"
OUT_DIR.mkdir(exist_ok=True)

# Canonical 4ct model results — NEVER use the v5 predecessor stale mirror at
# nmd/results/isoform_transitions/NMD_orf_model_v5/, which has different
# class ratios (4-CT scope drops to 22.4% NMD prevalence; v5 had 15.3%) and
# therefore different AUPRC headline numbers. See feedback memory
# "use the 4ct model, not the v5 predecessor mirror, for §5 figure work".
#
# Local clone of peter4244/NMD_orf_model_v5_4ct, results_4ct/ subdir.
# Files not present locally are SCP'd from the cluster on demand:
#   p.castaldi@explorer.northeastern.edu:/home/p.castaldi/cc/nmd_orf_model_v5_4ct/results_4ct/
MODEL_RESULTS = Path.home() / "claude_projects" / "NMD_orf_model_v5_4ct" / "results_4ct"

ATG_SIZES = [100, 500, 1000]
STOP_SIZES = [100, 500, 1000, 2000]
BEST_ATG = 500
BEST_STOP = 500
BEST_TAG = f"atg{BEST_ATG}_stop{BEST_STOP}"


def export_panel_d():
    """Best-model ROC curve — passes through (label, prob) per test transcript."""
    preds = pd.read_csv(MODEL_RESULTS / f"predictions_{BEST_TAG}.tsv", sep="\t")
    # Keep only what the ROC needs (smaller TSV downstream)
    out = OUT_DIR / "panelD_predictions.tsv"
    preds[["isoform_id", "label", "prob"]].to_csv(out, sep="\t", index=False)
    # And the headline metrics from the JSON
    with open(MODEL_RESULTS / f"metrics_{BEST_TAG}.json") as f:
        m = json.load(f)
    pd.DataFrame([{"tag": BEST_TAG, **m}]).to_csv(
        OUT_DIR / "panelD_metrics.tsv", sep="\t", index=False
    )
    n_test = int(m["n_test"])
    n_nmd = int(m["n_nmd"])
    print(f"  Panel D: wrote {out.name} ({len(preds):,} rows) + panelD_metrics.tsv")
    print(f"    Best model {BEST_TAG}: AUC={m['auc']:.4f}, AUPRC={m['auprc']:.4f}, "
          f"n_test={n_test:,} ({n_nmd:,} NMD)")


def export_panel_e():
    """Branch-level Shapley importance (Structural / Stop / ATG) on NMD test set.

    Mirrors orf_model_report fig2a_branch_importance: KernelSHAP with 500
    background, exact 2^3 coalitions over the three model branches.
    """
    df = pd.read_csv(MODEL_RESULTS / f"kernel_shap_branch_{BEST_TAG}.tsv", sep="\t")
    nmd = df[df["label"] > 0.5]

    rows = [
        {"branch": "Structural", "mean_abs_shap": nmd["shap_structural"].abs().mean()},
        {"branch": "Stop",       "mean_abs_shap": nmd["shap_stop"].abs().mean()},
        {"branch": "ATG",        "mean_abs_shap": nmd["shap_atg"].abs().mean()},
    ]
    out = pd.DataFrame(rows)
    total = out["mean_abs_shap"].sum()
    out["pct"] = 100 * out["mean_abs_shap"] / total
    out["n_nmd"] = len(nmd)

    path = OUT_DIR / "panelE_branch_importance.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"  Panel E: wrote {path.name} (n_nmd={len(nmd):,})")
    for _, r in out.iterrows():
        print(f"    {r['branch']:>12}: mean|SHAP|={r['mean_abs_shap']:.4f} ({r['pct']:.1f}%)")


# Friendly labels for the structural-branch channels (Panel F)
STRUCTURAL_LABELS = {
    "n_downstream_ejc": "Downstream EJC count",
    "is_sqanti_cds":    "Is TD2 CDS",
    "is_ref_cds":       "Is reference CDS",
    "frac_start":       "Fractional start position",
    "frac_stop":        "Fractional stop position",
}


def export_panel_f():
    """Per-feature DeepSHAP for the structural branch (averaged across 5 runs).

    Mirrors orf_model_report fig2b_structural_importance.
    """
    import glob
    files = sorted(glob.glob(str(MODEL_RESULTS / f"deepshap_summary_{BEST_TAG}_run*.tsv")))
    dfs = [pd.read_csv(f, sep="\t") for f in files]
    summary = (
        pd.concat(dfs)
        .groupby(["branch", "channel"], as_index=False)
        .agg({"mean_abs_shap_nmd": "mean"})
    )
    struct = (
        summary[summary["branch"] == "joint_structural"]
        .sort_values("mean_abs_shap_nmd", ascending=False)
        .copy()
    )
    struct["feature"] = struct["channel"].map(STRUCTURAL_LABELS).fillna(struct["channel"])
    struct = struct.rename(columns={"mean_abs_shap_nmd": "mean_abs_shap"})
    struct = struct[["feature", "channel", "mean_abs_shap"]]

    path = OUT_DIR / "panelF_structural_features.tsv"
    struct.to_csv(path, sep="\t", index=False)
    print(f"  Panel F: wrote {path.name} (n_runs={len(files)}, n_features={len(struct)})")
    for _, r in struct.iterrows():
        print(f"    {r['feature']:>28}: {r['mean_abs_shap']:.4f}")


def export_shap_profile(panel_letter, source_basename, out_name):
    """Per-position per-nucleotide |SHAP| for NMD and Control.

    Keeps the A/C/G/T channels separate (logomaker-ready), zoomed to ±20 nt
    around the codon (full ±250 nt available in the source).
    """
    df = pd.read_csv(MODEL_RESULTS / source_basename, sep="\t")
    nt = df[df["channel"].isin(["A", "C", "G", "T"])].copy()
    out = nt.sort_values(["relative_position", "channel"])[
        ["relative_position", "channel", "nmd_abs_shap", "ctrl_abs_shap"]
    ]
    path = OUT_DIR / out_name
    out.to_csv(path, sep="\t", index=False)
    print(f"  Panel {panel_letter}: wrote {path.name} "
          f"({len(out):,} rows, positions "
          f"{int(out.relative_position.min())}…{int(out.relative_position.max())})")
    return out


def export_panel_g():
    """ATG-window per-position |SHAP| (NMD vs Control, A+C+G+T sum)."""
    export_shap_profile("G", f"shap_profile_atg_joint_{BEST_TAG}.tsv",
                        "panelG_atg_profile.tsv")


def export_panel_h():
    """STOP-window per-position |SHAP| (NMD vs Control, A+C+G+T sum)."""
    export_shap_profile("H", f"shap_profile_stop_joint_{BEST_TAG}.tsv",
                        "panelH_stop_profile.tsv")


def main():
    print(f"Source: {MODEL_RESULTS}")
    export_panel_d()
    export_panel_e()
    export_panel_f()
    export_panel_g()
    export_panel_h()


if __name__ == "__main__":
    main()
