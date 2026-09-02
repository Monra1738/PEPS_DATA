"""Create per-crater APBT or CIRCLE summary CSVs from result JSONs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def crater_number(name: str) -> int:
    match = re.search(r"(\d+)", Path(str(name)).stem)
    if not match:
        raise ValueError(f"Cannot extract crater number from {name!r}")
    return int(match.group(1))


def _iter_literature_rows(value):
    if isinstance(value, list):
        return value
    if not isinstance(value, dict):
        return []
    for key in ("crater_database", "crater_candidates", "craters"):
        rows = value.get(key)
        if isinstance(rows, list):
            return rows
    return []


def load_locations(literature_path: Path) -> dict[int, tuple[object, object]]:
    if not literature_path.exists():
        return {}
    with literature_path.open() as fh:
        literature = json.load(fh)
    locations = {}
    for row in _iter_literature_rows(literature):
        crater_id = row.get("crater_id", row.get("id"))
        if crater_id is None:
            continue
        location = row.get("location") or {}
        locations[int(crater_id)] = (
            location.get("latitude"),
            location.get("longitude"),
        )
    return locations


def _finite_float(value, field: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"Non-finite value in {field}: {value!r}")
    return result


def apbt_run_values(run: dict) -> tuple[float, float, float, float]:
    diameter_max = _finite_float(
        run.get("Diameter_max_original", run["Diameter_max"]),
        "Diameter_max",
    )
    if run.get("Diameter_min") is not None:
        diameter_min = _finite_float(run["Diameter_min"], "Diameter_min")
    else:
        diameter_min = min(
            _finite_float(run["Diameter_ew"], "Diameter_ew"),
            _finite_float(run["Diameter_ns"], "Diameter_ns"),
        )
    diameter_mean = _finite_float(
        run.get("Diameter_mean", 0.5 * (diameter_min + diameter_max)),
        "Diameter_mean",
    )
    depth = run.get("Depth")
    if depth is None or float(depth) == 0.0:
        depth = run["Depth_med"]
    return diameter_mean, diameter_min, diameter_max, _finite_float(depth, "Depth")


def circle_run_values(run: dict) -> tuple[float, float]:
    diameter = _finite_float(run["Diameter_max"], "Diameter_max")
    depth = run.get("Depth")
    if depth is None or float(depth) == 0.0:
        depth = run["Depth_med"]
    return diameter, _finite_float(depth, "Depth")


def summarize_results(
    results: dict,
    *,
    method: str,
    locations: dict[int, tuple[object, object]] | None = None,
) -> list[dict]:
    method = method.upper()
    if method not in {"APBT", "CIRCLE"}:
        raise ValueError(f"Unsupported method: {method}")
    locations = locations or {}
    rows = []

    for crater_name in sorted(results, key=crater_number):
        runs = results[crater_name]
        if not isinstance(runs, list) or not runs:
            continue
        crater_id = crater_number(crater_name)
        latitude, longitude = locations.get(crater_id, (None, None))

        if method == "APBT":
            values = [apbt_run_values(run) for run in runs]
            mean_diameters, min_diameters, max_diameters, depths = map(
                np.asarray,
                zip(*values),
            )
            mean_diameter = float(np.mean(mean_diameters))
            min_diameter = float(np.min(min_diameters))
            max_diameter = float(np.max(max_diameters))
        else:
            values = [circle_run_values(run) for run in runs]
            diameters, depths = map(np.asarray, zip(*values))
            mean_diameter = float(np.mean(diameters))
            min_diameter = float(np.min(diameters))
            max_diameter = float(np.max(diameters))

        rows.append({
            "Crater": Path(crater_name).stem,
            "Lat": latitude,
            "Lon": longitude,
            "Mean Diameter": mean_diameter,
            "Max Diameter": max_diameter,
            "Min Diameter": min_diameter,
            "Mean Depth": float(np.mean(depths)),
            "Max Depth": float(np.max(depths)),
            "Min Depth": float(np.min(depths)),
        })
    return rows


FIELDNAMES = [
    "Crater",
    "Lat",
    "Lon",
    "Mean Diameter",
    "Max Diameter",
    "Min Diameter",
    "Mean Depth",
    "Max Depth",
    "Min Depth",
]


def write_summary(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def json_read(path: Path):
    with path.open() as fh:
        return json.load(fh)


def json_write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as fh:
        json.dump(value, fh, indent=2, default=float)
        fh.write("\n")
    temporary.replace(path)


def csv_write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def identity(crater: str, seed_points) -> tuple[str, tuple[int, ...]]:
    return Path(crater).stem.lower(), tuple(int(value) for value in seed_points)


def result_index(data: dict) -> tuple[dict, list[dict]]:
    index = {}
    duplicates = []
    for crater, runs in data.items():
        for run in runs:
            key = identity(crater, run["seed_points"])
            if key in index:
                duplicates.append({"crater": key[0], "seed_points": list(key[1])})
            index[key] = run
    return index, duplicates


def filter_results(
    data: dict,
    *,
    craters: set[str],
    identities: set[tuple[str, tuple[int, ...]]] | None = None,
) -> dict:
    selected = {}
    for crater, runs in data.items():
        crater_name = Path(crater).stem.lower()
        if crater_name not in craters:
            continue
        kept = [
            run
            for run in runs
            if identities is None or identity(crater_name, run["seed_points"]) in identities
        ]
        if kept:
            selected[crater_name] = kept
    return selected


def summary_rows(body: str, method: str, data: dict) -> list[dict]:
    return summarize_results(
        data,
        method=method,
        locations=load_locations(ROOT / "inputs" / body / "literature.json"),
    )


def case_metrics(case: dict, rows: list[dict]) -> dict:
    literature = case["literature_fn"]()
    values = []
    for row in rows:
        crater_id = crater_number(row["Crater"])
        reference = literature.get(crater_id)
        if not reference:
            continue
        literature_value = reference.get(case["literature_value_key"])
        if literature_value is None:
            continue
        measured = float(
            row["Mean Depth"] if case["metric"] == "depth" else row["Mean Diameter"]
        )
        reference_value = float(literature_value)
        error = measured - reference_value
        relative = abs(error) / abs(reference_value) * 100.0 if reference_value else math.nan
        literature_std = (
            reference.get(case["literature_std_key"])
            if case["literature_std_key"]
            else None
        )
        low = float(row["Min Depth"] if case["metric"] == "depth" else row["Min Diameter"])
        high = float(row["Max Depth"] if case["metric"] == "depth" else row["Max Diameter"])
        overlap = None
        if literature_std is not None and math.isfinite(float(literature_std)):
            lit_low = reference_value - float(literature_std)
            lit_high = reference_value + float(literature_std)
            overlap = not (high < lit_low or low > lit_high)
        values.append(
            {
                "crater": row["Crater"],
                "measured": measured,
                "literature": reference_value,
                "signed_error": error,
                "absolute_error": abs(error),
                "relative_error_pct": relative,
                "range_overlaps_1sigma": overlap,
            }
        )
    errors = np.asarray([value["signed_error"] for value in values], dtype=float)
    absolute = np.abs(errors)
    relative = np.asarray([value["relative_error_pct"] for value in values], dtype=float)
    overlap = [
        value["range_overlaps_1sigma"]
        for value in values
        if value["range_overlaps_1sigma"] is not None
    ]
    return {
        "n": len(values),
        "mean_absolute_error": float(np.mean(absolute)) if len(values) else None,
        "mape": float(np.mean(relative)) if len(values) else None,
        "rmse": float(np.sqrt(np.mean(errors**2))) if len(values) else None,
        "mean_signed_error_bias": float(np.mean(errors)) if len(values) else None,
        "strict_agreement_pct": float(np.mean(relative <= 10.0) * 100.0) if len(values) else None,
        "good_agreement_pct": float(np.mean(relative <= 25.0) * 100.0) if len(values) else None,
        "range_overlaps_1sigma_pct": float(np.mean(overlap) * 100.0) if overlap else None,
        "details": values,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--body", required=True)
    parser.add_argument("--method", choices=("APBT", "CIRCLE"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--literature", type=Path)
    return parser


def project_path(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    input_path = project_path(args.input)
    output_path = project_path(args.output)
    literature_path = project_path(
        args.literature or Path(f"inputs/{args.body}/literature.json")
    )
    with input_path.open() as fh:
        results = json.load(fh)
    rows = summarize_results(
        results,
        method=args.method,
        locations=load_locations(literature_path),
    )
    write_summary(rows, output_path)
    print(f"Saved {len(rows)} {args.body} {args.method} crater summaries to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
