#!/usr/bin/env python3
"""Build the verified publication dataset, comparisons, tables, and figures."""

from __future__ import annotations

import hashlib
import os
import shutil
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/apbt_publication_analysis_matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.batch_config import APBT_PATHS, CIRCLE_PATHS, project_path
from pipeline.summaries import (
    case_metrics,
    csv_write,
    filter_results,
    identity,
    json_read,
    json_write,
    result_index,
    summary_rows,
)
from pipeline.literature_validation import CASES
from pipeline.summaries import write_summary
from data.literature_bennu import return_lit_id_daly
from data.literature_didymos import return_lit_id_barnouin
from data.literature import compare_noguchi_hirata
from data.literature_bennu import compare_daly_bierhaus_deshapriya


BODIES = ("Bennu", "Ryugu", "Itokawa", "Didymos")
METHODS = ("APBT", "CIRCLE")
DEFAULT_OUTPUT = ROOT / "results" / "analysis"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def publication_sources() -> tuple[dict, dict]:
    """Return only the declared final result files.

    There are no historical-control fallbacks in the release pipeline.
    """
    gate = {
        "expansion_required": False,
        "policy": "use declared final results without rerunning measurements",
    }
    sources = {}
    for body in BODIES:
        sources[body] = {}
        for method in METHODS:
            defaults = APBT_PATHS[body] if method == "APBT" else CIRCLE_PATHS[body]
            path = project_path(defaults.output)
            if not path.is_file():
                raise FileNotFoundError(path)
            sources[body][method] = {
                "path": path,
                "provenance": "declared_final_release",
                "sha256": sha256(path),
            }
    return sources, gate


def expected_seed_count(body: str, seed_overrides: dict[str, Path] | None = None) -> int:
    seed_overrides = seed_overrides or {}
    seeds = json_read(seed_overrides.get(body, project_path(APBT_PATHS[body].seeds)))
    return sum(len(groups) for groups in seeds.values())


def build_publication_data(
    output: Path,
    sources: dict,
) -> dict:
    raw = {}
    summaries = {}
    manifest = {"policy": {}, "sources": {}}
    for body in BODIES:
        raw[body] = {}
        summaries[body] = {}
        manifest["sources"][body] = {}
        for method in METHODS:
            source = sources[body][method]
            data = json_read(source["path"])
            raw[body][method] = data
            rows = summary_rows(body, method, data)
            summaries[body][f"{method}_all_valid"] = rows
            summary_path = output / "results" / body / (
                "summary_table_per_crater.csv"
                if method == "APBT"
                else "circle_summary_table_per_crater.csv"
            )
            write_summary(rows, summary_path)
            valid_runs = len(result_index(data)[0])
            manifest["sources"][body][method] = {
                "path": str(source["path"].relative_to(ROOT)),
                "sha256": source["sha256"],
                "provenance": source["provenance"],
                "expected_seed_configurations": expected_seed_count(body),
                "valid_runs": valid_runs,
                "invalid_or_unavailable_runs": expected_seed_count(body) - valid_runs,
                "valid_craters": len(rows),
            }

        apbt_index = set(result_index(raw[body]["APBT"])[0])
        circle_index = set(result_index(raw[body]["CIRCLE"])[0])
        common_ids = apbt_index & circle_index
        for method in METHODS:
            paired_data = filter_results(
                raw[body][method],
                craters={Path(name).stem.lower() for name in raw[body][method]},
                identities=common_ids,
            )
            rows = summary_rows(body, method, paired_data)
            summaries[body][f"{method}_paired_common_seeds"] = rows
            write_summary(
                rows,
                output / "paired_summaries" / body.lower() / f"{method.lower()}_common_seeds.csv",
            )
        manifest["policy"][body] = {
            "standalone": "all available finite positive-depth results in the selected publication source",
            "direct_method_comparison": "identical APBT-CIRCLE valid seed identities",
            "paired_common_seed_count": len(common_ids),
        }

    json_write(output / "publication_dataset_manifest.json", manifest)
    return summaries


def build_literature_outputs(output: Path, summaries: dict, manifest: dict) -> list[dict]:
    metric_rows = []
    detail_rows = []
    for case in CASES:
        body = case["asteroid"]
        method = case["method"]
        datasets = {
            "standalone_all_valid": summaries[body][f"{method}_all_valid"],
            "paired_common_seeds": summaries[body][f"{method}_paired_common_seeds"],
        }
        for dataset, rows in datasets.items():
            metrics = case_metrics(case, rows)
            base = {
                "body": body,
                "method": method,
                "metric": case["metric"],
                "reference": case["reference"],
                "dataset": dataset,
                **{key: value for key, value in metrics.items() if key != "details"},
            }
            metric_rows.append(base)
            for detail in metrics["details"]:
                detail_rows.append({**base, **detail})
    csv_write(output / "comparisons" / "literature_validation_metrics.csv", metric_rows)
    csv_write(output / "comparisons" / "literature_validation_per_crater.csv", detail_rows)

    table_rows = []
    keys = []
    for row in metric_rows:
        if row["dataset"] != "standalone_all_valid":
            continue
        key = (row["body"], row["metric"], row["reference"])
        if key not in keys:
            keys.append(key)
    lookup = {
        (row["body"], row["method"], row["metric"], row["reference"]): row
        for row in metric_rows
        if row["dataset"] == "standalone_all_valid"
    }
    for body, metric, reference in keys:
        apbt = lookup.get((body, "APBT", metric, reference))
        circle = lookup.get((body, "CIRCLE", metric, reference))
        row = {"Target": body, "Quantity": metric.title(), "Reference": reference}
        for label, item in (("APBT", apbt), ("CIRCLE", circle)):
            source_info = manifest["sources"][body][label]
            row.update(
                {
                    f"{label} N": item["n"] if item else None,
                    f"{label} valid seeds": source_info["valid_runs"],
                    f"{label} invalid seeds": source_info["invalid_or_unavailable_runs"],
                    f"{label} 1sigma overlap pct": item["range_overlaps_1sigma_pct"] if item else None,
                    f"{label} strict pct": item["strict_agreement_pct"] if item else None,
                    f"{label} good pct": item["good_agreement_pct"] if item else None,
                    f"{label} MAE m": item["mean_absolute_error"] if item else None,
                    f"{label} MAPE pct": item["mape"] if item else None,
                    f"{label} RMSE m": item["rmse"] if item else None,
                    f"{label} bias m": item["mean_signed_error_bias"] if item else None,
                }
            )
        table_rows.append(row)
    csv_write(output / "manuscript_tables" / "table2_validation.csv", table_rows)
    return metric_rows


def save_literature_catalog_scatter(output: Path) -> None:
    image_dir = output / "figures" / "combined"
    image_dir.mkdir(parents=True, exist_ok=True)
    ryugu = pd.DataFrame.from_dict(compare_noguchi_hirata(), orient="index")
    bennu = pd.DataFrame.from_dict(compare_daly_bierhaus_deshapriya(), orient="index")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))

    ax = axes[0]
    ax.scatter(ryugu["hirata_diam"], ryugu["noguchi_diam"], alpha=0.72, s=42)
    limit = float(max(ryugu["hirata_diam"].max(), ryugu["noguchi_diam"].max())) * 1.05
    ax.plot([0, limit], [0, limit], "k--", linewidth=1.3)
    ax.set_xlabel("Crater diameters (m) [Hirata et al. (2020)]")
    ax.set_ylabel("Crater diameters (m) [Noguchi et al. (2021)]")
    ax.set_title("(a) Ryugu catalogue comparison")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.22)

    ax = axes[1]
    pair_specs = (
        ("daly_diam", "bierhaus_diam", "Daly vs Bierhaus"),
        ("daly_diam", "deshapriya_diam", "Daly vs Deshapriya"),
        ("bierhaus_diam", "deshapriya_diam", "Bierhaus vs Deshapriya"),
    )
    maximum = 0.0
    for x_col, y_col, label in pair_specs:
        subset = bennu[[x_col, y_col]].dropna()
        ax.scatter(subset[x_col], subset[y_col], alpha=0.66, s=38, label=label)
        if not subset.empty:
            maximum = max(maximum, float(subset.to_numpy().max()))
    limit = maximum * 1.05
    ax.plot([0, limit], [0, limit], "k--", linewidth=1.3, label="1:1 line")
    ax.set_xlabel("First catalogue diameter (m)")
    ax.set_ylabel("Second catalogue diameter (m)")
    ax.set_title("(b) Bennu pairwise catalogue comparisons")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        image_dir / "literature_diameter_catalogue_scatter.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def build_per_crater_manuscript_tables(output: Path, summaries: dict) -> None:
    def indexed(rows: list[dict]) -> dict[int, dict]:
        return {int(row["Crater"].rsplit("_", 1)[1]): row for row in rows}

    bennu_apbt = indexed(summaries["Bennu"]["APBT_all_valid"])
    bennu_circle = indexed(summaries["Bennu"]["CIRCLE_all_valid"])
    daly_ids = set(return_lit_id_daly())
    table3 = []
    for crater_id in sorted(set(bennu_apbt) - daly_ids):
        apbt = bennu_apbt[crater_id]
        circle = bennu_circle.get(crater_id)
        diameter = float(apbt["Mean Diameter"])
        depth = float(apbt["Mean Depth"])
        table3.append(
            {
                "Crater": crater_id,
                "Status": "not in Daly",
                "APBT D (m)": diameter,
                "APBT D min (m)": float(apbt["Min Diameter"]),
                "APBT D max (m)": float(apbt["Max Diameter"]),
                "APBT d (m)": depth,
                "APBT d/D": depth / diameter,
                "APBT d min (m)": float(apbt["Min Depth"]),
                "APBT d max (m)": float(apbt["Max Depth"]),
                "CIRCLE d (m)": float(circle["Mean Depth"]) if circle else None,
            }
        )
    csv_write(output / "manuscript_tables" / "table3_bennu_new_depths.csv", table3)

    didymos_apbt = indexed(summaries["Didymos"]["APBT_all_valid"])
    didymos_circle = indexed(summaries["Didymos"]["CIRCLE_all_valid"])
    barnouin_ids = set(return_lit_id_barnouin())
    table4 = []
    for crater_id in sorted(didymos_apbt):
        apbt = didymos_apbt[crater_id]
        circle = didymos_circle.get(crater_id)
        diameter = float(apbt["Mean Diameter"])
        depth = float(apbt["Mean Depth"])
        table4.append(
            {
                "Crater": crater_id,
                "Status": "diam. ref." if crater_id in barnouin_ids else "new candidate",
                "APBT D (m)": diameter,
                "APBT D min (m)": float(apbt["Min Diameter"]),
                "APBT D max (m)": float(apbt["Max Diameter"]),
                "APBT d (m)": depth,
                "APBT d/D": depth / diameter,
                "APBT d min (m)": float(apbt["Min Depth"]),
                "APBT d max (m)": float(apbt["Max Depth"]),
                "CIRCLE d (m)": float(circle["Mean Depth"]) if circle else None,
            }
        )
    csv_write(output / "manuscript_tables" / "table4_didymos_morphometry.csv", table4)


def run_existing_apbt_cross_asteroid(output: Path) -> None:
    from pipeline import scaling as cross

    cross.ROOT = output
    cross.OUT_DIR = output / "results" / "cross_asteroid"
    cross.CSV_DIR = cross.OUT_DIR / "csv"
    cross.IMG_DIR = cross.OUT_DIR / "images"
    cross.main()


def run_combined_validation_figures(output: Path) -> None:
    from pipeline import figures as combined

    combined.ROOT = output
    combined.OUTPUT_DIR = output / "figures" / "combined"
    combined.main()


def write_analysis_summary(output: Path, sources: dict, gate: dict, metrics: list[dict]) -> None:
    standalone = [row for row in metrics if row["dataset"] == "standalone_all_valid"]
    lines = [
        "# Verified presented-paper analysis",
        "",
        f"- Bennu/Ryugu full rerun triggered: {gate['expansion_required']}",
        f"- Literature-validation cases: {len(standalone)}",
        "",
        "## Standalone literature validation",
        "",
        "| Body | Method | Metric | Reference | N | MAE (m) | MAPE (%) |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for row in standalone:
        lines.append(
            f"| {row['body']} | {row['method']} | {row['metric']} | {row['reference']} | "
            f"{row['n']} | {row['mean_absolute_error']:.3f} | {row['mape']:.3f} |"
        )
    lines.extend(["", "## Source provenance", ""])
    for body in BODIES:
        for method in METHODS:
            source = sources[body][method]
            lines.append(
                f"- {body} {method}: `{source['path'].relative_to(ROOT)}` "
                f"({source['provenance']})."
            )
    lines.append("")
    (output / "analysis_report.md").write_text("\n".join(lines))


def main() -> None:
    """Rebuild all derived CSV and PNG files."""

    output = DEFAULT_OUTPUT
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)

    print("[1/3] Reading final APBT and CIRCLE measurements")
    sources, gate = publication_sources()
    summaries = build_publication_data(output, sources)
    manifest = json_read(output / "publication_dataset_manifest.json")

    print("[2/3] Comparing measurements with literature data")
    metric_rows = build_literature_outputs(output, summaries, manifest)

    print("[3/3] Creating CSV tables and PNG figures")
    build_per_crater_manuscript_tables(output, summaries)
    save_literature_catalog_scatter(output)
    run_existing_apbt_cross_asteroid(output)
    run_combined_validation_figures(output)
    write_analysis_summary(output, sources, gate, metric_rows)
    print(f"Created analysis outputs in {output}")


if __name__ == "__main__":
    main()
