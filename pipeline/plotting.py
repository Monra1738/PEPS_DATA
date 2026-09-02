"""Small shared publication-plot style."""

from __future__ import annotations

import matplotlib.pyplot as plt


PLOT_FONT_SIZES = {
    "title": 18,
    "label": 16,
    "tick": 13,
    "legend": 13,
    "annotation": 11,
}


def configure_publication_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": PLOT_FONT_SIZES["tick"],
            "axes.titlesize": PLOT_FONT_SIZES["title"],
            "axes.labelsize": PLOT_FONT_SIZES["label"],
            "xtick.labelsize": PLOT_FONT_SIZES["tick"],
            "ytick.labelsize": PLOT_FONT_SIZES["tick"],
            "legend.fontsize": PLOT_FONT_SIZES["legend"],
            "figure.titlesize": PLOT_FONT_SIZES["title"],
        }
    )

