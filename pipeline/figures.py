import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/apbt_matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from data.literature import return_lit_id_hirata, return_lit_id_noguchi
from data.literature_bennu import (
    return_lit_id_bierhaus,
    return_lit_id_daly,
    return_lit_id_deshapriya,
)
from data.literature_didymos import return_lit_id_barnouin
from data.literature_itokawa import return_lit_id_naru_hirata
from pipeline.plotting import PLOT_FONT_SIZES, configure_publication_plot_style
from utils.ids import extract_crater_id

OUTPUT_DIR = ROOT / "results" / "cross_asteroid" / "images"


def method_label(method: str) -> str:
    return "APBT v2" if method == "APBT" else "CIRCLE"


def load_summary(asteroid: str, method: str) -> pd.DataFrame:
    prefix = "circle_" if method == "CIRCLE" else ""
    path = ROOT / "results" / asteroid / f"{prefix}summary_table_per_crater.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path).set_index("Crater")


def build_literature_comparison(
    *,
    asteroid: str,
    method: str,
    metric: str,
    literature_data: dict[int, dict],
    literature_key: str,
) -> pd.DataFrame:
    summary = load_summary(asteroid, method)
    value_col = "Mean Depth" if metric == "depth" else "Mean Diameter"
    rows = []
    for crater_name, row in summary.iterrows():
        crater_id = extract_crater_id(crater_name)
        literature = literature_data.get(crater_id)
        if not literature:
            continue

        literature_value = literature.get(literature_key)
        measured_value = row.get(value_col)
        if pd.isna(literature_value) or pd.isna(measured_value):
            continue

        rows.append({
            "crater": crater_name,
            "literature": float(literature_value),
            "measured": float(measured_value),
        })

    return pd.DataFrame(rows)


def build_didymos_depth_comparison() -> pd.DataFrame:
    apbt = load_summary("Didymos", "APBT")
    circle = load_summary("Didymos", "CIRCLE")
    df = pd.DataFrame({
        "literature": apbt["Mean Depth"],
        "measured": circle["Mean Depth"],
    }).dropna()
    df["crater"] = df.index
    return df.reset_index(drop=True)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.05,
        f"({label})",
        transform=ax.transAxes,
        fontsize=PLOT_FONT_SIZES["title"],
        fontweight="bold",
        va="bottom",
        ha="left",
        clip_on=False,
    )


def plot_scatter_panel(
    ax,
    df: pd.DataFrame,
    *,
    panel_label: str,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    if df.empty:
        raise ValueError(f"No data for panel {panel_label}: {title}")

    x = df["literature"].to_numpy(dtype=float)
    y = df["measured"].to_numpy(dtype=float)
    ax.scatter(
        x,
        y,
        color="#4c72b0",
        alpha=0.72,
        edgecolors="white",
        linewidths=0.6,
        s=56,
        label="Craters",
    )

    lim_min = 0.0
    lim_max = float(max(x.max(), y.max())) * 1.1
    if lim_max == 0:
        lim_max = 1.0

    ax.plot(
        [lim_min, lim_max],
        [lim_min, lim_max],
        color="black",
        linestyle="--",
        linewidth=1.5,
        alpha=0.65,
        label="1:1 line",
    )

    if len(df) >= 2:
        slope, intercept_value = np.polyfit(x, y, 1)
        x_plot = np.array([lim_min, lim_max])
        intercept_text = (
            f"- {abs(intercept_value):.2f}"
            if intercept_value < 0
            else f"+ {intercept_value:.2f}"
        )
        ax.plot(
            x_plot,
            slope * x_plot + intercept_value,
            color="red",
            linestyle="-",
            linewidth=2,
            label=f"fit: y = {slope:.2f}x {intercept_text}",
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_xlabel(xlabel, fontsize=PLOT_FONT_SIZES["label"])
    ax.set_ylabel(ylabel, fontsize=PLOT_FONT_SIZES["label"])
    ax.set_title(title, fontsize=PLOT_FONT_SIZES["title"], fontweight="bold")
    ax.tick_params(axis="both", labelsize=PLOT_FONT_SIZES["tick"])
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(frameon=True, loc="upper left", fontsize=PLOT_FONT_SIZES["legend"])
    add_panel_label(ax, panel_label)


def comparison_panel(
    *,
    asteroid: str,
    method: str,
    metric: str,
    reference_label: str,
    literature_data: dict[int, dict],
    literature_key: str,
    panel_label: str,
) -> dict:
    df = build_literature_comparison(
        asteroid=asteroid,
        method=method,
        metric=metric,
        literature_data=literature_data,
        literature_key=literature_key,
    )
    display_method = method_label(method)
    return {
        "df": df,
        "panel_label": panel_label,
        "title": f"{display_method} vs {reference_label}",
        "xlabel": f"{reference_label} {metric} (m)",
        "ylabel": f"{display_method} {metric} (m)",
    }


def didymos_depth_panel(panel_label: str) -> dict:
    return {
        "df": build_didymos_depth_comparison(),
        "panel_label": panel_label,
        "title": "CIRCLE vs APBT v2 depth",
        "xlabel": "APBT v2 depth (m)",
        "ylabel": "CIRCLE depth (m)",
    }


def save_combined_figure(
    panels: list[dict],
    *,
    filename: str,
    ncols: int,
    figsize: tuple[float, float],
) -> None:
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    axes_array = np.atleast_1d(axes).ravel()

    for ax, panel in zip(axes_array, panels):
        plot_scatter_panel(ax, **panel)

    for ax in axes_array[len(panels):]:
        ax.set_visible(False)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / filename
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


def main() -> None:
    configure_publication_plot_style()

    noguchi = return_lit_id_noguchi()
    hirata_ryugu = return_lit_id_hirata()
    daly = return_lit_id_daly()
    bierhaus = return_lit_id_bierhaus()
    deshapriya = return_lit_id_deshapriya()
    hirata_itokawa = return_lit_id_naru_hirata()
    barnouin = return_lit_id_barnouin()

    save_combined_figure(
        [
            comparison_panel(
                asteroid="Ryugu",
                method="APBT",
                metric="diameter",
                reference_label="Noguchi et al. (2021)",
                literature_data=noguchi,
                literature_key="noguchi_diameter_mean",
                panel_label="a",
            ),
            comparison_panel(
                asteroid="Ryugu",
                method="CIRCLE",
                metric="diameter",
                reference_label="Noguchi et al. (2021)",
                literature_data=noguchi,
                literature_key="noguchi_diameter_mean",
                panel_label="b",
            ),
            comparison_panel(
                asteroid="Ryugu",
                method="APBT",
                metric="diameter",
                reference_label="Hirata et al. (2020)",
                literature_data=hirata_ryugu,
                literature_key="hirata_diameter_m",
                panel_label="c",
            ),
            comparison_panel(
                asteroid="Ryugu",
                method="CIRCLE",
                metric="diameter",
                reference_label="Hirata et al. (2020)",
                literature_data=hirata_ryugu,
                literature_key="hirata_diameter_m",
                panel_label="d",
            ),
        ],
        filename="ryugu_diameter_validation.png",
        ncols=2,
        figsize=(12.6, 12.2),
    )

    save_combined_figure(
        [
            comparison_panel(
                asteroid="Ryugu",
                method="APBT",
                metric="depth",
                reference_label="Noguchi et al. (2021)",
                literature_data=noguchi,
                literature_key="noguchi_depth_mean",
                panel_label="a",
            ),
            comparison_panel(
                asteroid="Ryugu",
                method="CIRCLE",
                metric="depth",
                reference_label="Noguchi et al. (2021)",
                literature_data=noguchi,
                literature_key="noguchi_depth_mean",
                panel_label="b",
            ),
        ],
        filename="ryugu_depth_validation.png",
        ncols=2,
        figsize=(12.6, 6.0),
    )

    save_combined_figure(
        [
            comparison_panel(
                asteroid="Bennu",
                method="APBT",
                metric="diameter",
                reference_label="Daly et al. (2020)",
                literature_data=daly,
                literature_key="daly_diameter_mean",
                panel_label="a",
            ),
            comparison_panel(
                asteroid="Bennu",
                method="CIRCLE",
                metric="diameter",
                reference_label="Daly et al. (2020)",
                literature_data=daly,
                literature_key="daly_diameter_mean",
                panel_label="b",
            ),
            comparison_panel(
                asteroid="Bennu",
                method="APBT",
                metric="depth",
                reference_label="Daly et al. (2020)",
                literature_data=daly,
                literature_key="daly_depth_mean",
                panel_label="c",
            ),
            comparison_panel(
                asteroid="Bennu",
                method="CIRCLE",
                metric="depth",
                reference_label="Daly et al. (2020)",
                literature_data=daly,
                literature_key="daly_depth_mean",
                panel_label="d",
            ),
        ],
        filename="bennu_daly_validation.png",
        ncols=2,
        figsize=(12.6, 12.2),
    )

    save_combined_figure(
        [
            comparison_panel(
                asteroid="Bennu",
                method="APBT",
                metric="diameter",
                reference_label="Bierhaus et al. (2022)",
                literature_data=bierhaus,
                literature_key="bierhaus_diameter_mean",
                panel_label="a",
            ),
            comparison_panel(
                asteroid="Bennu",
                method="CIRCLE",
                metric="diameter",
                reference_label="Bierhaus et al. (2022)",
                literature_data=bierhaus,
                literature_key="bierhaus_diameter_mean",
                panel_label="b",
            ),
            comparison_panel(
                asteroid="Bennu",
                method="APBT",
                metric="diameter",
                reference_label="Deshapriya et al. (2021)",
                literature_data=deshapriya,
                literature_key="deshapriya_diameter_m",
                panel_label="c",
            ),
            comparison_panel(
                asteroid="Bennu",
                method="CIRCLE",
                metric="diameter",
                reference_label="Deshapriya et al. (2021)",
                literature_data=deshapriya,
                literature_key="deshapriya_diameter_m",
                panel_label="d",
            ),
        ],
        filename="bennu_bierhaus_deshapriya_diameter_validation.png",
        ncols=2,
        figsize=(12.6, 12.2),
    )

    save_combined_figure(
        [
            comparison_panel(
                asteroid="Itokawa",
                method="APBT",
                metric="diameter",
                reference_label="Hirata et al. (2009)",
                literature_data=hirata_itokawa,
                literature_key="naru_hirata_diameter_m",
                panel_label="a",
            ),
            comparison_panel(
                asteroid="Itokawa",
                method="CIRCLE",
                metric="diameter",
                reference_label="Hirata et al. (2009)",
                literature_data=hirata_itokawa,
                literature_key="naru_hirata_diameter_m",
                panel_label="b",
            ),
            comparison_panel(
                asteroid="Itokawa",
                method="APBT",
                metric="depth",
                reference_label="Hirata et al. (2009)",
                literature_data=hirata_itokawa,
                literature_key="naru_hirata_depth_m",
                panel_label="c",
            ),
            comparison_panel(
                asteroid="Itokawa",
                method="CIRCLE",
                metric="depth",
                reference_label="Hirata et al. (2009)",
                literature_data=hirata_itokawa,
                literature_key="naru_hirata_depth_m",
                panel_label="d",
            ),
        ],
        filename="itokawa_validation.png",
        ncols=2,
        figsize=(12.6, 12.2),
    )

    save_combined_figure(
        [
            comparison_panel(
                asteroid="Didymos",
                method="APBT",
                metric="diameter",
                reference_label="Barnouin et al. (2024)",
                literature_data=barnouin,
                literature_key="Barnouin_diameter_m",
                panel_label="a",
            ),
            comparison_panel(
                asteroid="Didymos",
                method="CIRCLE",
                metric="diameter",
                reference_label="Barnouin et al. (2024)",
                literature_data=barnouin,
                literature_key="Barnouin_diameter_m",
                panel_label="b",
            ),
            didymos_depth_panel("c"),
        ],
        filename="didymos_validation.png",
        ncols=3,
        figsize=(18.6, 6.0),
    )


if __name__ == "__main__":
    main()
