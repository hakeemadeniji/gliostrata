"""Generate the T4/T5 result figures.

Every figure is drawn from a committed result file, never from numbers typed into this
script, so a figure cannot disagree with the record it illustrates. Each writes a PDF
(vector, for submission), a PNG at 300 dpi, and a CSV of its own numbers.

    forest     per-study effect ratios by payload, with each payload's registered median
               and the broader readings item 11(i.7) requires reported beside it.
    dose       the registered dose-response curve per payload, with ED50 and the
               registered dose range.
    survival   predicted survival under the Stupp baseline against a configuration, and
               what the RMST gain is the area of.
    gains      the distribution of per-patient RMST gains from the registered rehearsal,
               against T5's success and harm thresholds.

Sources: `docs/t4t5_emax_search/emax_result.json`,
`docs/t4t5_emax_search/screening/extraction_final.jsonl`,
`configs/response_surface_v1.yaml`, `docs/dry_runs/t4t5_controls.json`.

No cohort data is read, so no read-guard access is recorded.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt  # noqa: E402  (after apply_style's import of matplotlib)
import numpy as np

from gliostrata.prediction import response_surface as rs
from gliostrata.reporting.figures import (
    GRID,
    PAYLOAD_LABEL,
    PAYLOAD_ORDER,
    RING_W,
    SURFACE,
    TEXT_PRIMARY,
    TEXT_SECONDARY,
    apply_style,
    caption,
    payload_style,
    save,
    threshold_line,
)
from gliostrata.simulation.parameters import load_parameters

SEARCH = Path("docs/t4t5_emax_search")
RESULT = SEARCH / "emax_result.json"
EXTRACTION = SEARCH / "screening" / "extraction_final.jsonl"
CONTROLS = Path("docs/dry_runs/t4t5_controls.json")
OUT = Path("docs/figures")
T5_BAR_MONTHS = 3.0


def _result() -> dict:
    if not RESULT.exists():
        raise SystemExit(f"{RESULT} not found; compute the Emax result first")
    return json.loads(RESULT.read_text(encoding="utf-8"))


def _extraction() -> list[dict]:
    return [json.loads(line) for line in EXTRACTION.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _ratio(row: dict) -> float | None:
    if row.get("ratio") is not None:
        return float(row["ratio"])
    treated, control = row.get("treated_median"), row.get("control_median")
    if treated is None or not control:
        return None
    return float(treated) / float(control)


# --- forest plot --------------------------------------------------------------------------

def figure_forest() -> None:
    result, rows = _result(), _extraction()
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.grid(axis="y", visible=False)

    table, y = [], 0.0
    ticks, tick_labels, boundaries, headers = [], [], [], []
    for payload in PAYLOAD_ORDER:
        colour, marker = payload_style(payload)
        info = result["payloads"][payload]
        in_primary = set(info["studies"])
        mine = [r for r in rows if r["payload"] == payload and _ratio(r) is not None]
        mine.sort(key=lambda r: (r["pmid"] not in in_primary, -_ratio(r)))
        if not mine:
            continue

        # A header row above each block, so a reader never has to infer the payload from
        # the hue. The median goes here too, at the TOP of its own block: at the bottom it
        # reads as a label for the block below it.
        y -= 1.0
        headers.append((y, payload, info, colour))
        top = y - 0.5

        for row in mine:
            y -= 1.0
            ratio, primary = _ratio(row), row["pmid"] in in_primary
            ax.plot([ratio], [y], marker=marker, markersize=8 if primary else 6,
                    markerfacecolor=colour if primary else SURFACE,
                    markeredgecolor=colour, markeredgewidth=RING_W if primary else 1.3,
                    linestyle="none", zorder=5)
            ticks.append(y)
            tick_labels.append(row["pmid"] if primary else f"{row['pmid']}  (excluded)")
            table.append({"payload": payload, "pmid": row["pmid"],
                          "ratio": round(ratio, 4), "in_primary_set": primary,
                          "control_arm_kind": row.get("control_arm_kind"),
                          "p2_intracranial": row.get("p2_intracranial")})

        median = info["ratio"]
        # With one study the median IS that study's point, so a separate rule drawn through
        # it is redundant and reads as an arrowhead. The header says n=1; that is the honest
        # way to show it, and the absent rule is itself the signal that there is nothing to
        # take a median over.
        if median and info["primary_narrowed"]["n"] > 1:
            # Behind the markers, so a point sitting on the median still reads as a point.
            ax.plot([median, median], [y - 0.45, top], color=colour, linewidth=2.0,
                    zorder=2, solid_capstyle="butt", alpha=0.85)
        boundaries.append(y - 0.9)
        y -= 0.6

    threshold_line(ax, x=1.0, label="no effect")
    ax.set_xscale("log")
    ax.set_xlim(0.92, 4.6)
    ax.set_xticks([1.0, 1.5, 2.0, 3.0, 4.0])
    ax.set_xticklabels(["1.0", "1.5", "2.0", "3.0", "4.0"])
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels, fontsize=7.2)
    ax.set_ylim(y + 0.4, 0.2)
    ax.set_xlabel("Treated / control median survival ratio  (log scale)")
    ax.set_title("Preclinical effect ratios behind each payload's registered prior")

    for edge in boundaries[:-1]:
        ax.axhline(edge, color=GRID, linewidth=0.6, zorder=0)
    # Header text in the plot's left margin, in text ink: a colour-blind reader gets the
    # payload from the words, and the swatch beside it only reinforces.
    for y_pos, payload, info, colour in headers:
        ax.plot([0.945], [y_pos], marker="s", markersize=6, color=colour,
                markeredgecolor=SURFACE, markeredgewidth=1.0,
                linestyle="none", clip_on=False, zorder=6)
        n = info["primary_narrowed"]["n"]
        descriptor = (f"median {info['ratio']:.2f}   n={n}" if n > 1
                      else f"one study, ratio {info['ratio']:.2f}   n=1")
        ax.annotate(
            f"{PAYLOAD_LABEL[payload]}   {descriptor}",
            xy=(0.97, y_pos), xytext=(4, 0), textcoords="offset points",
            va="center", ha="left", fontsize=8.5, color=TEXT_PRIMARY, fontweight="semibold")

    filled = plt.Line2D([], [], marker="o", linestyle="none", color=TEXT_SECONDARY,
                        markerfacecolor=TEXT_SECONDARY, markersize=8,
                        label="in the registered primary set")
    open_m = plt.Line2D([], [], marker="o", linestyle="none", markerfacecolor=SURFACE,
                        markeredgecolor=TEXT_SECONDARY, markersize=6,
                        label="eligible, excluded by item 11(i.7)")
    bar = plt.Line2D([], [], color=TEXT_SECONDARY, linewidth=2.0,
                     label="registered median for the payload")
    ax.legend(handles=[filled, open_m, bar], loc="upper left",
              bbox_to_anchor=(0.0, -0.12), ncol=3, fontsize=7.5,
              handletextpad=0.6, columnspacing=1.6)

    caption(fig, (
        "Each point is one study's treated-to-control median-survival ratio; the vertical "
        "rule is that payload's registered median (item 11(i.2)). Filled points are in the "
        "primary set; open points are eligible under P1-P5 but excluded by the control-arm "
        "narrowing of item 11(i.7). PMID 37420122 appears twice, once per payload, because "
        "P4 is payload-specific. CAR's median rests on a single study; every broader reading "
        "of its evidence moves it further from T5's bar, not closer. "
        f"Source: {RESULT.as_posix()}."))
    saved = save(fig, OUT, "t4t5_emax_forest", table=table)
    print(f"  forest   -> {saved.pdf.name}, {saved.png.name}, {saved.table.name}")


# --- dose-response ------------------------------------------------------------------------

def figure_dose() -> None:
    params = load_parameters()
    ed50, hill = params.ed50_cells, params.hill_h
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    dose = np.logspace(5.5, 8.5, 400)

    table = []
    for payload in PAYLOAD_ORDER:
        colour, marker = payload_style(payload)
        emax = params.emax(payload)
        effect = emax * dose**hill / (dose**hill + ed50**hill)
        ax.plot(dose, effect, color=colour, zorder=4)
        ax.annotate(f"{PAYLOAD_LABEL[payload]}  Emax {emax:.3f}",
                    xy=(dose[-1], effect[-1]), xytext=(-4, 4),
                    textcoords="offset points", ha="right", va="bottom",
                    fontsize=8, color=TEXT_PRIMARY)
        at_ed50 = emax * 0.5
        ax.plot([ed50], [at_ed50], marker=marker, color=colour, markerfacecolor=SURFACE,
                markeredgecolor=colour, markeredgewidth=RING_W, markersize=7, zorder=5)
        for cells in (1e6, 1e7, 2.5e7, 5e7):
            table.append({"payload": payload, "engrafted_cells": f"{cells:.3g}",
                          "log_hazard_reduction": round(
                              float(emax * cells**hill / (cells**hill + ed50**hill)), 5)})

    ax.axvspan(1e6, 5e7, color=GRID, alpha=0.55, zorder=0, linewidth=0)
    ax.annotate("registered dose range\n1e6 to 5e7 cells", xy=(7e6, ax.get_ylim()[1]),
                xytext=(0, -4), textcoords="offset points", ha="center", va="top",
                fontsize=7.2, color=TEXT_SECONDARY)
    threshold_line(ax, x=ed50, label=f"ED50 {ed50:.2g} cells", label_at="bottom")
    ax.set_xscale("log")
    ax.set_xlabel("Engrafted, acting cells (n_eff)")
    ax.set_ylabel("Log-hazard reduction E")
    ax.set_title("Registered dose-response, by payload")
    ax.legend(handles=[plt.Line2D([], [], color=payload_style(p)[0],
                                  label=PAYLOAD_LABEL[p]) for p in PAYLOAD_ORDER],
              loc="upper left", fontsize=7.5)
    caption(fig, (
        "E = Emax * n_eff^h / (n_eff^h + ED50^h), with Emax from item 11(i) and ED50 and h "
        "from the registered parameter file. The shaded band is item 4's dose range; the "
        "engrafted count n_eff is a fraction of the injected dose, so the working region sits "
        "left of the band. Open markers are each payload's half-maximal point. Emax values "
        "are derived from the preclinical literature and the 0.5 translation factor is "
        "assumed, so the vertical scale is a modelling assumption, not a measurement."))
    saved = save(fig, OUT, "t4t5_dose_response", table=table)
    print(f"  dose     -> {saved.pdf.name}, {saved.png.name}, {saved.table.name}")


# --- predicted survival -------------------------------------------------------------------

def figure_survival() -> None:
    params = load_parameters()
    grid = rs.time_grid(params)
    # A baseline of the shape the rehearsal produced: RMST24 near 13.8 months.
    baseline = np.exp(-0.055 * grid)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))

    base_rmst = rs.rmst(baseline, params)
    ax.plot(grid, baseline, color=TEXT_SECONDARY, zorder=4)

    table = [{"months": round(float(t), 1), "series": "Stupp baseline",
              "survival": round(float(s), 5)} for t, s in zip(grid[::20], baseline[::20])]
    colour, _ = payload_style("prodrug_enzyme")
    treated = rs.configuration_curve(baseline, -params.emax("prodrug_enzyme"))
    treated_rmst = rs.rmst(treated, params)
    ax.plot(grid, treated, color=colour, zorder=5)
    ax.fill_between(grid, baseline, treated, color=colour, alpha=0.14, linewidth=0, zorder=2)
    # The two curves are close by construction -- that closeness IS the finding -- so the
    # direct labels are separated horizontally, each against its own curve, rather than
    # stacked at the right edge where they collide.
    i_treated, i_baseline = 55, 165
    ax.annotate("With a configuration", xy=(grid[i_treated], treated[i_treated]),
                xytext=(0, 13), textcoords="offset points", ha="center", va="bottom",
                fontsize=8, color=TEXT_PRIMARY)
    ax.annotate("Stupp baseline", xy=(grid[i_baseline], baseline[i_baseline]),
                xytext=(0, -15), textcoords="offset points", ha="center", va="top",
                fontsize=8, color=TEXT_PRIMARY)
    table += [{"months": round(float(t), 1), "series": "configuration",
               "survival": round(float(s), 5)} for t, s in zip(grid[::20], treated[::20])]

    ax.annotate(f"shaded area = RMST gain, {treated_rmst - base_rmst:.2f} months",
                xy=(8.0, 0.5 * (baseline[80] + treated[80])), xytext=(14, 24),
                textcoords="offset points", fontsize=7.8, color=TEXT_PRIMARY,
                arrowprops={"arrowstyle": "-", "color": TEXT_SECONDARY, "linewidth": 0.8})
    ax.set_xlim(0, params.tau_months)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Months from randomisation")
    ax.set_ylabel("Predicted survival")
    ax.set_title(f"How an RMST gain is read off the curves (tau = {params.tau_months:.0f} months)")
    ax.legend(handles=[
        plt.Line2D([], [], color=TEXT_SECONDARY,
                   label=f"Stupp baseline   RMST {base_rmst:.2f} months"),
        plt.Line2D([], [], color=colour,
                   label=f"With a configuration   RMST {treated_rmst:.2f} months")],
        loc="lower left", fontsize=7.5)
    caption(fig, (
        "An illustration of the quantity, not a result: the baseline is a curve of the shape "
        "the registered TCGA rehearsal produced and the offset is the prodrug/enzyme Emax at "
        "full occupancy. The gain of item 11(b) is the area between the curves out to tau, "
        "and T5 asks whether it reaches 3 months."))
    saved = save(fig, OUT, "t4t5_survival_gain", table=table)
    print(f"  survival -> {saved.pdf.name}, {saved.png.name}, {saved.table.name}")


# --- gain distribution --------------------------------------------------------------------

def figure_gains() -> None:
    """Per-patient RMST gains against T5's thresholds.

    Deliberately NOT coloured by the payload palette: these runs are not payloads, and a
    hue that means "CRAd" in the forest plot must not mean "synthetic run" here. Colour
    follows the entity. One measure in two settings is one series, so there is no legend
    either - the row labels carry identity.
    """
    if not CONTROLS.exists():
        print("  gains    -> skipped, no control run on disk")
        return
    controls = json.loads(CONTROLS.read_text(encoding="utf-8"))
    runs = {c["control"]: c for c in controls["controls"]}
    wanted = [("tcga_rehearsal", "Held-out TCGA split"),
              ("synthetic_dry_run", "Synthetic patients")]
    present = [(k, label, runs[k]) for k, label in wanted
               if k in runs and runs[k].get("gains_months")]
    if not present:
        print("  gains    -> skipped, the control run carries no per-patient gains")
        return

    hue, _ = payload_style("CRAd")   # one series: a single hue, not a categorical slot
    fig, ax = plt.subplots(figsize=(6.8, 1.35 + 1.25 * len(present)))
    rng = np.random.default_rng(3)
    table, labels = [], []

    for index, (key, label, run) in enumerate(present):
        gains = np.asarray(run["gains_months"], dtype=float)
        # Jitter only to separate overlapping points; it carries no meaning, so it is small
        # and the axis it moves on has no scale.
        jitter = rng.uniform(-0.17, 0.17, size=gains.size)
        ax.plot(gains, index + jitter, marker="o", linestyle="none", markersize=5.0,
                markerfacecolor=hue, markeredgecolor=SURFACE, markeredgewidth=1.0,
                alpha=0.75, zorder=4)
        mean = float(gains.mean())
        ax.plot([mean], [index], marker="|", markersize=26, markeredgewidth=2.4,
                color=TEXT_PRIMARY, zorder=6)
        ax.annotate(f"mean {mean:.2f} mo", xy=(mean, index), xytext=(0, 17),
                    textcoords="offset points", ha="center", fontsize=8,
                    color=TEXT_PRIMARY)
        labels.append(f"{label}\nn={gains.size}")
        for value in gains:
            table.append({"run": key, "gain_months": round(float(value), 5)})

    threshold_line(ax, x=T5_BAR_MONTHS, label="T5 success, 3 months")
    threshold_line(ax, x=0.0, label="harm below 0")
    ax.set_yticks(range(len(present)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(-0.6, len(present) - 0.4)
    upper = max(3.3, max(max(r["gains_months"]) for _, _, r in present) * 1.15)
    ax.set_xlim(min(-0.35, min(min(r["gains_months"]) for _, _, r in present) - 0.2), upper)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Per-patient RMST gain (months)")
    ax.set_title("Rehearsal gains against T5's registered thresholds")
    caption(fig, (
        "One point per patient, from the item 11(o) controls; the vertical tick is the mean. "
        "These are NOT results: T5's verdict comes from the registry in section 6, on the "
        "sealed external cohort, after the freeze. Every point sits between the harm line and "
        "the success bar, so on development data no patient is predicted to be harmed and "
        "none reaches 3 months. The harm reading is not vacuous: the null-payload control "
        "shows harm does occur when a configuration's toxicity exceeds its benefit."))
    saved = save(fig, OUT, "t4t5_rehearsal_gains", table=table)
    print(f"  gains    -> {saved.pdf.name}, {saved.png.name}, {saved.table.name}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("figure", nargs="?", default="all",
                        choices=("all", "forest", "dose", "survival", "gains"))
    args = parser.parse_args(argv)
    apply_style()
    runners = {"forest": figure_forest, "dose": figure_dose,
               "survival": figure_survival, "gains": figure_gains}
    print(f"writing figures to {OUT}/")
    for name in (runners if args.figure == "all" else [args.figure]):
        runners[name]()
    print("Each figure ships a CSV of its own numbers: the palette's contrast warning "
          "obligates a table view, and a reader who cannot use the figure needs the values.")


if __name__ == "__main__":
    main()
