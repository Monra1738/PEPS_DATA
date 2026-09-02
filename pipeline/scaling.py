from __future__ import annotations

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

from pipeline.plotting import PLOT_FONT_SIZES, configure_publication_plot_style

configure_publication_plot_style()

ASTEROIDS = ("Bennu", "Ryugu", "Itokawa", "Didymos")
OUT_DIR = ROOT / "results" / "cross_asteroid"
CSV_DIR = OUT_DIR / "csv"
IMG_DIR = OUT_DIR / "images"

ASTEROID_COLORS = {
    "Bennu": "#1f77b4",
    "Ryugu": "#2ca02c",
    "Itokawa": "#d62728",
    "Didymos": "#9467bd",
}


def crater_id_from_name(name: str) -> int:
    return int(str(name).split("_")[1])


def load_summary_table(asteroid: str) -> pd.DataFrame:
    path = ROOT / "results" / asteroid / "summary_table_per_crater.csv"
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    df = df.rename(columns={
        "Crater": "crater",
        "Lat": "lat",
        "Lon": "lon",
        "Mean Diameter": "diameter_m",
        "Max Diameter": "diameter_max_m",
        "Min Diameter": "diameter_min_m",
        "Mean Depth": "depth_m",
        "Max Depth": "depth_max_m",
        "Min Depth": "depth_min_m",
    })
    df["asteroid"] = asteroid
    df["crater_id"] = df["crater"].map(crater_id_from_name)

    numeric_cols = [
        "diameter_m",
        "diameter_min_m",
        "diameter_max_m",
        "depth_m",
        "depth_min_m",
        "depth_max_m",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["diameter_range_m"] = df["diameter_max_m"] - df["diameter_min_m"]
    df["depth_range_m"] = df["depth_max_m"] - df["depth_min_m"]
    df["d_to_D"] = df["depth_m"] / df["diameter_m"]
    df["d_to_D_min"] = df["depth_min_m"] / df["diameter_max_m"]
    df["d_to_D_max"] = df["depth_max_m"] / df["diameter_min_m"]
    df["diameter_range_pct"] = 100.0 * df["diameter_range_m"] / df["diameter_m"]
    df["depth_range_pct"] = 100.0 * df["depth_range_m"] / df["depth_m"]

    return df[[
        "asteroid",
        "crater_id",
        "crater",
        "lat",
        "lon",
        "diameter_m",
        "diameter_min_m",
        "diameter_max_m",
        "diameter_range_m",
        "diameter_range_pct",
        "depth_m",
        "depth_min_m",
        "depth_max_m",
        "depth_range_m",
        "depth_range_pct",
        "d_to_D",
        "d_to_D_min",
        "d_to_D_max",
    ]]


def load_all_measurements() -> pd.DataFrame:
    frames = []
    for asteroid in ASTEROIDS:
        frames.append(load_summary_table(asteroid))
    df = pd.concat(frames, ignore_index=True)
    return df.sort_values(["asteroid", "crater_id"]).reset_index(drop=True)


def iqr(series: pd.Series) -> float:
    return float(series.quantile(0.75) - series.quantile(0.25))


def grouped_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for asteroid, group in df.groupby("asteroid", sort=True):
        rows.append({
            "asteroid": asteroid,
            "n_craters": len(group),
            "diameter_min_m": group["diameter_m"].min(),
            "diameter_median_m": group["diameter_m"].median(),
            "diameter_max_m": group["diameter_m"].max(),
            "depth_min_m": group["depth_m"].min(),
            "depth_median_m": group["depth_m"].median(),
            "depth_max_m": group["depth_m"].max(),
            "d_to_D_mean": group["d_to_D"].mean(),
            "d_to_D_median": group["d_to_D"].median(),
            "d_to_D_iqr": iqr(group["d_to_D"]),
            "d_to_D_min": group["d_to_D"].min(),
            "d_to_D_max": group["d_to_D"].max(),
            "diameter_range_pct_median": group["diameter_range_pct"].median(),
            "depth_range_pct_median": group["depth_range_pct"].median(),
        })
    return pd.DataFrame(rows)


def fit_loglog(group: pd.DataFrame, label: str) -> dict[str, float | str | int]:
    work = group[
        (group["diameter_m"] > 0)
        & (group["depth_m"] > 0)
        & group["diameter_m"].notna()
        & group["depth_m"].notna()
    ].copy()

    row: dict[str, float | str | int] = {"fit": label, "n_craters": len(work)}
    if len(work) < 3:
        row.update({
            "coefficient_a": np.nan,
            "exponent_b": np.nan,
            "r2_loglog": np.nan,
        })
        return row

    log_diameter = np.log10(work["diameter_m"].to_numpy(dtype=float))
    log_depth = np.log10(work["depth_m"].to_numpy(dtype=float))
    exponent_b, intercept = np.polyfit(log_diameter, log_depth, 1)
    predicted_log = exponent_b * log_diameter + intercept
    ss_res = float(np.sum((log_depth - predicted_log) ** 2))
    ss_tot = float(np.sum((log_depth - log_depth.mean()) ** 2))
    r2_log = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan

    row.update({
        "coefficient_a": 10 ** intercept,
        "exponent_b": exponent_b,
        "r2_loglog": r2_log,
    })
    return row


def scaling_fit_summary(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([fit_loglog(df, "APBT_all_asteroids")])


def add_global_scaling_residuals(df: pd.DataFrame, fits: pd.DataFrame) -> pd.DataFrame:
    fit = fits[fits["fit"] == "APBT_all_asteroids"].iloc[0]
    coefficient = float(fit["coefficient_a"])
    exponent = float(fit["exponent_b"])
    out = df.copy()
    out["global_scaling_predicted_depth_m"] = coefficient * out["diameter_m"] ** exponent
    out["depth_ratio_to_global_scaling"] = (
        out["depth_m"] / out["global_scaling_predicted_depth_m"]
    )
    out["log10_depth_residual"] = (
        np.log10(out["depth_m"]) - np.log10(out["global_scaling_predicted_depth_m"])
    )
    return out


def residual_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for asteroid, group in df.groupby("asteroid", sort=True):
        rows.append({
            "asteroid": asteroid,
            "n_craters": len(group),
            "median_depth_ratio_to_global_scaling": group[
                "depth_ratio_to_global_scaling"
            ].median(),
            "mean_depth_ratio_to_global_scaling": group[
                "depth_ratio_to_global_scaling"
            ].mean(),
            "median_log10_depth_residual": group["log10_depth_residual"].median(),
            "iqr_log10_depth_residual": iqr(group["log10_depth_residual"]),
            "pct_deeper_than_global_scaling": (
                group["depth_ratio_to_global_scaling"] > 1.0
            ).mean() * 100.0,
        })
    return pd.DataFrame(rows)


def common_diameter_range(df: pd.DataFrame) -> tuple[float, float]:
    mins = df.groupby("asteroid")["diameter_m"].min()
    maxs = df.groupby("asteroid")["diameter_m"].max()
    return float(mins.max()), float(maxs.min())


def common_range_summary(df: pd.DataFrame) -> pd.DataFrame:
    low, high = common_diameter_range(df)
    common = df[(df["diameter_m"] >= low) & (df["diameter_m"] <= high)].copy()
    rows = []
    for asteroid, group in common.groupby("asteroid", sort=True):
        rows.append({
            "asteroid": asteroid,
            "common_diameter_min_m": low,
            "common_diameter_max_m": high,
            "n_craters": len(group),
            "diameter_median_m": group["diameter_m"].median(),
            "depth_median_m": group["depth_m"].median(),
            "d_to_D_median": group["d_to_D"].median(),
            "depth_ratio_to_global_scaling_median": group[
                "depth_ratio_to_global_scaling"
            ].median(),
            "log10_depth_residual_median": group["log10_depth_residual"].median(),
        })
    return pd.DataFrame(rows)


def table5_summary(
    summary: pd.DataFrame,
    residuals: pd.DataFrame,
    common_range: pd.DataFrame,
) -> pd.DataFrame:
    summary_by_body = summary.set_index("asteroid")
    residuals_by_body = residuals.set_index("asteroid")
    common_by_body = common_range.set_index("asteroid")
    rows = []
    for asteroid in ASTEROIDS:
        base = summary_by_body.loc[asteroid]
        residual = residuals_by_body.loc[asteroid]
        common = common_by_body.loc[asteroid]
        rows.append(
            {
                "Asteroid": asteroid,
                "n": int(base["n_craters"]),
                "Median D (m)": float(base["diameter_median_m"]),
                "Median d (m)": float(base["depth_median_m"]),
                "Median d/D": float(base["d_to_D_median"]),
                "Median d/pred.": float(
                    residual["median_depth_ratio_to_global_scaling"]
                ),
                "Common-range d/pred.": float(
                    common["depth_ratio_to_global_scaling_median"]
                ),
            }
        )
    return pd.DataFrame(rows)


def save_residual_boxplot(df: pd.DataFrame) -> None:
    labels = []
    data = []
    colors = []
    for asteroid in ASTEROIDS:
        group = df[df["asteroid"] == asteroid]
        if group.empty:
            continue
        labels.append(asteroid)
        data.append(group["depth_ratio_to_global_scaling"].dropna().to_numpy(dtype=float))
        colors.append(ASTEROID_COLORS[asteroid])

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    box = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.4)
    ax.set_ylabel("Depth / global APBT scaling prediction")
    ax.set_title("Size-normalized APBT depth residuals")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(
        IMG_DIR / "apbt_cross_asteroid_global_scaling_residual_boxplot.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_common_range_residual_plot(common_summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.2))
    x = np.arange(len(common_summary))
    colors = [ASTEROID_COLORS[a] for a in common_summary["asteroid"]]
    ax.bar(
        x,
        common_summary["depth_ratio_to_global_scaling_median"],
        color=colors,
        alpha=0.72,
    )
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.4)
    ax.set_xticks(x)
    ax.set_xticklabels(common_summary["asteroid"])
    ax.set_ylabel("Median depth / global scaling prediction")
    low = common_summary["common_diameter_min_m"].iloc[0]
    high = common_summary["common_diameter_max_m"].iloc[0]
    ax.set_title(f"APBT residuals in common diameter range ({low:.1f}-{high:.1f} m)")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(
        IMG_DIR / "apbt_cross_asteroid_common_range_residuals.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_scatter_loglog(df: pd.DataFrame, fits: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))

    for asteroid in ASTEROIDS:
        group = df[df["asteroid"] == asteroid]
        if group.empty:
            continue
        ax.scatter(
            group["diameter_m"],
            group["depth_m"],
            label=asteroid,
            c=ASTEROID_COLORS[asteroid],
            alpha=0.72,
            edgecolors="white",
            linewidths=0.4,
            s=42,
        )

    x_min = df["diameter_m"].min() * 0.85
    x_max = df["diameter_m"].max() * 1.15
    x_line = np.logspace(np.log10(x_min), np.log10(x_max), 200)
    fit = fits[fits["fit"] == "APBT_all_asteroids"]
    if not fit.empty and pd.notna(fit.iloc[0]["coefficient_a"]):
        row = fit.iloc[0]
        y_line = row["coefficient_a"] * x_line ** row["exponent_b"]
        ax.plot(
            x_line,
            y_line,
            linestyle="-",
            color="black",
            linewidth=2,
            label=f"APBT scaling b={row['exponent_b']:.2f}",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Diameter D (m)")
    ax.set_ylabel("Depth d (m)")
    ax.set_title("APBT cross-asteroid crater depth-diameter scaling")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=PLOT_FONT_SIZES["legend"], ncols=2)
    fig.tight_layout()
    fig.savefig(
        IMG_DIR / "apbt_cross_asteroid_depth_diameter_loglog.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    df = load_all_measurements()
    summary = grouped_summary(df)
    fits = scaling_fit_summary(df)
    df = add_global_scaling_residuals(df, fits)
    residuals = residual_summary(df)
    common_range = common_range_summary(df)
    table5 = table5_summary(summary, residuals, common_range)

    df.to_csv(CSV_DIR / "apbt_cross_asteroid_per_crater_metrics.csv", index=False)
    summary.to_csv(CSV_DIR / "apbt_cross_asteroid_d_to_D_summary.csv", index=False)
    fits.to_csv(CSV_DIR / "apbt_cross_asteroid_scaling_fits.csv", index=False)
    residuals.to_csv(CSV_DIR / "apbt_cross_asteroid_global_scaling_residual_summary.csv", index=False)
    common_range.to_csv(CSV_DIR / "apbt_cross_asteroid_common_diameter_range_summary.csv", index=False)
    table_dir = ROOT / "manuscript_tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    table5.to_csv(table_dir / "table5_cross_asteroid_summary.csv", index=False)

    save_scatter_loglog(df, fits)
    save_residual_boxplot(df)
    save_common_range_residual_plot(common_range)

    print("Saved cross-asteroid CSV outputs:")
    for path in sorted(CSV_DIR.glob("*.csv")):
        print(f"- {path.relative_to(ROOT)}")
    print("\nSaved cross-asteroid figures:")
    for path in sorted(IMG_DIR.glob("*.png")):
        print(f"- {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
