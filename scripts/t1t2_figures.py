"""Result figures for T1 and T2.

Both tests have been evaluated, so unlike the T4/T5 figures these show results rather than
rehearsals, and each carries its registered verdict. Every number is read from a committed
result file; none is typed into this script.

    t1-deconvolution   estimated against realised state fraction, one panel per cellular
                       state, over the 100 registered mixtures. Registered PASS.
    t2-external        each variant's discrimination and calibration on CGGA against their
                       registered thresholds. Registered FAIL on both rows.

Sources: `data/results/t1/neftel_mixture_deconvolution.json` and its mixture table,
`data/results/cgga/cox_external_validation.json`, `data/results/cgga/brier18_external.json`.

**No cohort data is read.** `data/results/cgga/` holds the summary statistics already
published in the pre-registration, not the sealed cohort under `data/raw/cgga` or
`data/processed/cgga`, which this script never opens.

The four source files were untracked when these figures were first written, although the
pre-registration cites all four by path, including in two section 8 amendment rows. A
registration that names a file has to ship it, or the citation cannot be checked, so they are
now committed as narrow exceptions to the `data/results/**` ignore. Everything else under
`data/results/` stays ignored.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd

from gliostrata.reporting.figures import (
    GRID,
    RING_W,
    SERIES,
    SURFACE,
    TEXT_SECONDARY,
    apply_style,
    caption,
    save,
    threshold_line,
)

T1_RESULT = Path("data/results/t1/neftel_mixture_deconvolution.json")
T1_MIXTURES = Path("data/results/t1/neftel_mixture_deconvolution_mixtures.parquet")
T2_COX = Path("data/results/cgga/cox_external_validation.json")
T2_BRIER = Path("data/results/cgga/brier18_external.json")
OUT = Path("docs/figures")

STATES = ("NPC", "OPC", "AC", "MES")
T1_THRESHOLD = 0.75
T2_CINDEX_THRESHOLD = 0.65
T2_BRIER_THRESHOLD = 0.25
VARIANTS = (
    ("variant_1_baseline_committed", "Variant 1\nage, sex, KPS, arm"),
    ("variant_2_transferable", "Variant 2\nage, sex, arm"),
    ("variant_3_molecular", "Variant 3\n+ IDH, MGMT"),
)


# --- T1: deconvolution ---------------------------------------------------------------------

def figure_t1() -> None:
    result = json.loads(T1_RESULT.read_text(encoding="utf-8"))
    mixtures = pd.read_parquet(T1_MIXTURES)
    primary = result["primary"]
    per_state = primary["per_state_r"]
    boot = primary["bootstrap"]["per_state"]

    # Small multiples of ONE measure, so one hue across all four panels. Four categorical
    # hues would put four slots on screen together, which the palette does not support on
    # the all-pairs list - and the state is already named in each panel's title.
    hue = SERIES[0]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.6), sharex=True, sharey=True)
    table = []

    # One scale across all four panels, computed once over every state. Small multiples are
    # only comparable on a common scale, and with shared axes a per-panel limit would be
    # silently overwritten by whichever panel was drawn last.
    hi = 1.08 * max(
        max(mixtures[f"realised_{s}"].max(), mixtures[f"estimate_{s}"].max())
        for s in STATES)
    lo = 0.0

    for ax, state in zip(axes.ravel(), STATES):
        realised = mixtures[f"realised_{state}"].to_numpy(dtype=float)
        estimate = mixtures[f"estimate_{state}"].to_numpy(dtype=float)
        ax.plot([lo, hi], [lo, hi], color=TEXT_SECONDARY, linewidth=1.0,
                linestyle=(0, (4, 3)), zorder=2)
        ax.plot(realised, estimate, marker="o", linestyle="none", markersize=4.4,
                markerfacecolor=hue, markeredgecolor=SURFACE, markeredgewidth=0.8,
                alpha=0.8, zorder=4)
        r = per_state[state]
        ax.set_title(f"{state}    r = {r:.3f}", fontsize=9)
        ax.annotate(f"95% CI {boot[state]['lower']:.3f} to {boot[state]['upper']:.3f}",
                    xy=(0.04, 0.95), xycoords="axes fraction", va="top", ha="left",
                    fontsize=7.4, color=TEXT_SECONDARY)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        for value_r, value_e in zip(realised, estimate):
            table.append({"state": state, "realised_fraction": round(float(value_r), 6),
                          "estimated_fraction": round(float(value_e), 6)})

    for ax in axes[-1]:
        ax.set_xlabel("Realised fraction")
    for ax in axes[:, 0]:
        ax.set_ylabel("Estimated fraction")

    mean_r = primary["mean_r"]
    mean_ci = primary["bootstrap"]["mean_r"]
    fig.suptitle(
        f"T1 state deconvolution: mean r = {mean_r:.4f}  "
        f"(95% CI {mean_ci['lower']:.4f} to {mean_ci['upper']:.4f})   PASS vs "
        f"{T1_THRESHOLD}",
        x=0.0, ha="left", fontsize=10.5, fontweight="semibold", y=1.0)
    caption(fig, (
        f"One point per mixture, {len(mixtures)} registered pseudo-bulk mixtures from the "
        "Neftel 2019 atlas (pre-registration item 10, spec in docs/t1_deconvolution_spec.md). "
        "The dashed line is identity, where a perfect estimate would sit. The registered "
        f"metric is the mean of the four per-state Pearson r values against a threshold of "
        f"{T1_THRESHOLD}, decided on the 95% CI lower bound; it PASSES. This row reads a "
        "simulated reference, not a patient cohort, so it tests the estimator rather than "
        "its transfer to real tissue."))
    fig.tight_layout()
    saved = save(fig, OUT, "t1_deconvolution", table=table)
    print(f"  t1   -> {saved.pdf.name}, {saved.png.name}, {saved.table.name}")


# --- T2: external validation ----------------------------------------------------------------

def figure_t2() -> None:
    cox = json.loads(T2_COX.read_text(encoding="utf-8"))
    brier = json.loads(T2_BRIER.read_text(encoding="utf-8"))

    # Two measures on different scales and in opposite directions, so two panels. Never two
    # y-scales on one axis.
    fig, (ax_c, ax_b) = plt.subplots(1, 2, figsize=(8.4, 3.5), sharey=True)
    internal_hue, external_hue = SERIES[0], SERIES[1]
    table = []

    for index, (key, _label) in enumerate(VARIANTS):
        y = len(VARIANTS) - 1 - index
        entry = cox[key]
        internal = entry.get("internal_cindex")
        if internal is not None:
            ax_c.plot([internal], [y], marker="o", markersize=7, markerfacecolor=SURFACE,
                      markeredgecolor=internal_hue, markeredgewidth=RING_W,
                      linestyle="none", zorder=5)
        if entry.get("external_evaluable"):
            point = entry["external_cindex"]
            lower, upper = entry["external_ci_lower"], entry["external_ci_upper"]
            ax_c.plot([lower, upper], [y, y], color=external_hue, linewidth=2.0, zorder=4,
                      solid_capstyle="butt")
            ax_c.plot([point], [y], marker="s", markersize=7, color=external_hue,
                      markeredgecolor=SURFACE, markeredgewidth=1.2, linestyle="none",
                      zorder=6)
            table.append({"variant": key, "metric": "external C-index",
                          "point": round(point, 6), "ci_lower": round(lower, 6),
                          "ci_upper": round(upper, 6),
                          "threshold": T2_CINDEX_THRESHOLD, "clears": False})
        else:
            ax_c.annotate("NOT EVALUABLE\nCGGA has no KPS", xy=(0.52, y), xytext=(0, 0),
                          textcoords="offset points", fontsize=7.4, va="center",
                          ha="left", color=TEXT_SECONDARY)
        if internal is not None:
            table.append({"variant": key, "metric": "internal C-index",
                          "point": round(internal, 6), "ci_lower": None, "ci_upper": None,
                          "threshold": None, "clears": None})

        b_entry = brier.get(key, {})
        if b_entry.get("external_evaluable"):
            point = b_entry["brier_18mo"]
            lower, upper = b_entry["ci_lower"], b_entry["ci_upper"]
            ax_b.plot([lower, upper], [y, y], color=external_hue, linewidth=2.0, zorder=4,
                      solid_capstyle="butt")
            ax_b.plot([point], [y], marker="s", markersize=7, color=external_hue,
                      markeredgecolor=SURFACE, markeredgewidth=1.2, linestyle="none",
                      zorder=6)
            table.append({"variant": key, "metric": "Brier @ 18 months",
                          "point": round(point, 6), "ci_lower": round(lower, 6),
                          "ci_upper": round(upper, 6),
                          "threshold": T2_BRIER_THRESHOLD, "clears": False})
        else:
            ax_b.annotate("NOT EVALUABLE", xy=(0.205, y), fontsize=7.4, va="center",
                          ha="left", color=TEXT_SECONDARY)

    threshold_line(ax_c, x=T2_CINDEX_THRESHOLD, label=f"needs >= {T2_CINDEX_THRESHOLD}")
    threshold_line(ax_b, x=T2_BRIER_THRESHOLD, label=f"needs <= {T2_BRIER_THRESHOLD}")
    ax_c.set_xlim(0.45, 0.72)
    ax_b.set_xlim(0.20, 0.30)
    for ax, title in ((ax_c, "Discrimination (C-index)"),
                      (ax_b, "Calibration (Brier at 18 months)")):
        ax.set_yticks(range(len(VARIANTS)))
        ax.set_yticklabels([label for _, label in reversed(VARIANTS)], fontsize=7.6)
        ax.set_ylim(-0.6, len(VARIANTS) - 0.4)
        ax.grid(axis="y", visible=False)
        ax.set_title(title, fontsize=9.5)
        ax.axhline(0.5, color=GRID, linewidth=0.6, zorder=0)
        ax.axhline(1.5, color=GRID, linewidth=0.6, zorder=0)

    ax_c.legend(handles=[
        plt.Line2D([], [], marker="o", linestyle="none", markerfacecolor=SURFACE,
                   markeredgecolor=internal_hue, markeredgewidth=RING_W, markersize=7,
                   label="Internal, TCGA (cross-validated)"),
        plt.Line2D([], [], marker="s", color=external_hue, markersize=7,
                   label="External, CGGA (95% CI)")],
        loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=2, fontsize=7.5)

    fig.suptitle("T2 external validation on CGGA: both registered rows FAIL",
                 x=0.0, ha="left", fontsize=10.5, fontweight="semibold", y=1.02)
    caption(fig, (
        "Registered thresholds are a C-index of at least 0.65 on the CI lower bound and a "
        "Brier score at 18 months of at most 0.25 on the CI upper bound; the verdict is "
        "decided by variant 2, the transferable model. Both rows FAIL. Variant 1, the "
        "committed baseline, is NOT EVALUABLE externally because CGGA records no KPS, which "
        "is a property of the cohort and not a result. The internal values are "
        "cross-validated on TCGA and are shown only to make the drop on transfer visible; "
        "they are not a second test. Variant 3 carries a harmonisation caveat: CGGA's IDH "
        "and MGMT calling is not necessarily identical to TCGA's."))
    fig.tight_layout()
    saved = save(fig, OUT, "t2_external_validation", table=table)
    print(f"  t2   -> {saved.pdf.name}, {saved.png.name}, {saved.table.name}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("figure", nargs="?", default="all",
                        choices=("all", "t1", "t2"))
    args = parser.parse_args(argv)
    apply_style()
    print(f"writing figures to {OUT}/")
    if args.figure in ("all", "t1"):
        figure_t1()
    if args.figure in ("all", "t2"):
        figure_t2()


if __name__ == "__main__":
    main()
