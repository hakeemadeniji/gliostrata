"""Publication figure style for GlioStrata results.

One validated visual system, applied by every figure in the project, so a colour means
the same thing in every panel and no script picks its own.

**The palette is computed, not chosen.** The three categorical slots below are the first
three of a documented eight-hue theme, and they were checked with that theme's validator on
the all-pairs list, which is the right list here because these figures put all three series
on screen together (forest plots, overlaid curves):

    slots 1-3, light surface, --pairs all
      lightness band        PASS, all inside L 0.43-0.78
      chroma floor          PASS, all >= 0.1
      CVD separation        PASS, worst pair dE 12.6 (deutan), 10.8 (tritan); target >= 8
      normal-vision floor   PASS, worst pair dE 23.0; floor >= 15
      contrast vs surface   WARN, orange 2.88 and aqua 2.54 are below 3:1

That WARN carries an obligation rather than a choice: a palette below 3:1 must ship
**visible direct labels or a table view**. Both are provided here - every figure direct-labels
its series, and `write_table` emits the figure's own numbers as CSV beside it. Each series also
gets a distinct marker, so identity never rests on colour alone in print, under
colour-blindness, or in a photocopy.

Rules this module enforces for its callers, each of which is a known way a chart goes wrong:

* categorical hues are assigned **in fixed order and never cycled** - a fourth payload
  would not get a generated hue, it would need a facet;
* grid and axis rules are **solid hairlines**, one shade off the surface. Dashing is
  reserved for a real threshold, where it means something;
* text wears text inks, never a series colour;
* **no dual axis, ever** - two scales become two panels;
* marks are thin: 1.8pt lines, >= 6pt markers, a surface-coloured ring where marks overlap;
* values are direct-labelled **selectively**, never one number per point.

Figures are written as PDF (vector, for submission) and PNG at 300 dpi (for preview).
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # no display in CI or on a headless run
import matplotlib.pyplot as plt  # noqa: E402

# -- the design system's parameters, the only values a brand swap would touch ---------------

SURFACE = "#ffffff"
TEXT_PRIMARY = "#191919"
TEXT_SECONDARY = "#666666"
GRID = "#e6e6e6"          # one shade off the surface, solid
AXIS = "#cccccc"

# Categorical slots 1-3 of the documented theme, in its fixed order.
SERIES = ("#257de0", "#f47230", "#00b894")       # blue, orange, aqua
# Secondary encoding, so identity survives loss of colour.
MARKERS = ("o", "s", "^")

PAYLOAD_ORDER = ("CRAd", "prodrug_enzyme", "CAR")
PAYLOAD_LABEL = {"CRAd": "CRAd", "prodrug_enzyme": "Prodrug/enzyme", "CAR": "CAR"}

LINE_W = 1.8
MARKER_SIZE = 7.0
RING_W = 1.6               # surface-coloured ring on overlapping marks
HAIRLINE = 0.6


def payload_style(payload: str) -> tuple[str, str]:
    """The colour and marker for a payload, by its fixed position. Never cycled.

    A payload outside the registered three is an error rather than a fourth generated hue:
    the validated palette covers three slots on the all-pairs list, and a fourth would have
    to be a facet.
    """
    try:
        index = PAYLOAD_ORDER.index(payload)
    except ValueError:
        raise ValueError(
            f"{payload!r} is not one of the registered payloads {PAYLOAD_ORDER}. The "
            "validated palette covers exactly these three on the all-pairs list; a fourth "
            "series needs a facet, not a new hue.") from None
    return SERIES[index], MARKERS[index]


def apply_style() -> None:
    """Set the rcParams every figure in the project shares."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.labelsize": 9,
        "axes.labelcolor": TEXT_SECONDARY,
        "axes.edgecolor": AXIS,
        "axes.linewidth": HAIRLINE,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": HAIRLINE,
        "grid.linestyle": "-",          # solid: dashing is reserved for a real threshold
        "xtick.color": TEXT_SECONDARY,
        "ytick.color": TEXT_SECONDARY,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": 8,
        "lines.linewidth": LINE_W,
        "lines.markersize": MARKER_SIZE,
        "lines.solid_capstyle": "round",
        "text.color": TEXT_PRIMARY,
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,             # embed as TrueType, not Type 3: journals require it
        "ps.fonttype": 42,
    })


def threshold_line(ax, *, y: float | None = None, x: float | None = None,
                   label: str | None = None, color: str = TEXT_SECONDARY,
                   label_at: str = "top") -> None:
    """A registered threshold or null value. Dashed, because here dashing MEANS something."""
    if (y is None) == (x is None):
        raise ValueError("threshold_line takes exactly one of y or x")
    if y is not None:
        ax.axhline(y, color=color, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    else:
        ax.axvline(x, color=color, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    if label:
        if y is not None:
            ax.annotate(label, xy=(1.0, y), xycoords=("axes fraction", "data"),
                        xytext=(-2, 3), textcoords="offset points",
                        ha="right", va="bottom", fontsize=7.5, color=TEXT_SECONDARY)
        else:
            # `label_at` exists because a top-anchored threshold label collides with the
            # direct labels at the end of rising curves. Moving the label is the fix; moving
            # the curve's label away from its curve would not be.
            at_top = label_at == "top"
            ax.annotate(label, xy=(x, 1.0 if at_top else 0.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(3, -2 if at_top else 3), textcoords="offset points",
                        ha="left", va="top" if at_top else "bottom",
                        fontsize=7.5, color=TEXT_SECONDARY)


def caption(fig, text: str) -> None:
    """A figure's provenance and caveats, in text ink, below the plot."""
    fig.text(0.0, -0.02, text, ha="left", va="top", fontsize=7.2,
             color=TEXT_SECONDARY, wrap=True)


@dataclass(frozen=True)
class Saved:
    pdf: Path
    png: Path
    table: Path | None


def save(fig, out_dir: Path, stem: str, *, table: Sequence[dict[str, Any]] | None = None,
         table_fields: Sequence[str] | None = None) -> Saved:
    """Write the figure as vector and raster, and its numbers as a table.

    The table is not optional decoration. The palette's contrast WARN obligates a table
    view or visible labels, and a reader who cannot use the figure must still be able to
    read the values.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf, png = out_dir / f"{stem}.pdf", out_dir / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png)
    plt.close(fig)

    table_path = None
    if table:
        table_path = out_dir / f"{stem}.csv"
        fields = list(table_fields or table[0].keys())
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(table)
    return Saved(pdf=pdf, png=png, table=table_path)
