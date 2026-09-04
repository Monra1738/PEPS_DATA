"""Simple staged reviewer workflow for APBT v2 and CIRCLE."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import tempfile
from importlib import import_module
from pathlib import Path

from pipeline.literature_validation import CASES
from pipeline.summaries import (
    case_metrics,
    csv_write,
    filter_results,
    json_read,
    result_index,
    summary_rows,
    write_summary,
)
from scripts.batch_config import APBT_PATHS, CIRCLE_PATHS, ROOT, project_path


BODIES = ("Bennu", "Ryugu", "Itokawa", "Didymos")
METHODS = ("APBT", "CIRCLE")
RUNS_ROOT = ROOT / "runs"
RUN_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class WorkflowError(RuntimeError):
    """A short, actionable workflow error."""


def _read_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, default=float)
        handle.write("\n")
    os.replace(temporary, path)


def _portable(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _stored_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _run_dir(name: str) -> Path:
    if not RUN_NAME.fullmatch(name):
        raise WorkflowError(
            "Run names may contain only letters, numbers, dots, dashes, and underscores."
        )
    return RUNS_ROOT / name


def _manifest(run_dir: Path, name: str) -> dict:
    path = run_dir / "run.json"
    if path.exists():
        return _read_json(path)
    return {"version": 1, "name": name, "results": {}}


def _save_manifest(run_dir: Path, manifest: dict) -> None:
    _write_json(run_dir / "run.json", manifest)


def _body(value: str) -> str:
    if value.lower() == "all":
        return "all"
    lookup = {body.lower(): body for body in BODIES}
    try:
        return lookup[value.lower()]
    except KeyError as exc:
        raise argparse.ArgumentTypeError(
            f"unknown body {value!r}; choose {', '.join(BODIES)}, or all"
        ) from exc


def _selected_bodies(values: list[str] | None, manifest: dict) -> list[str]:
    if values:
        return list(BODIES) if "all" in values else list(dict.fromkeys(values))
    saved = manifest.get("selection", {}).get("bodies")
    return list(saved) if saved else list(BODIES)


def _selected_methods(value: str | None, manifest: dict) -> list[str]:
    if value:
        return list(METHODS) if value == "both" else [value.upper()]
    saved = manifest.get("selection", {}).get("methods")
    return list(saved) if saved else list(METHODS)


def _result_path(manifest: dict, body: str, method: str) -> Path:
    try:
        value = manifest["results"][body][method]["path"]
    except KeyError as exc:
        raise WorkflowError(
            f"No {method} measurements for {body}. Run the measure stage first."
        ) from exc
    path = _stored_path(value)
    if not path.is_file():
        raise WorkflowError(f"Measurement file is missing: {path}")
    return path


def _remove_appledouble(root: Path) -> None:
    if not root.exists():
        return
    for path in root.rglob("._*"):
        if path.is_file():
            path.unlink()


def _invalidate_derived(run_dir: Path, manifest: dict) -> None:
    comparison = run_dir / "comparison.json"
    if comparison.exists():
        comparison.unlink()
    outputs = run_dir / "outputs"
    if outputs.exists():
        shutil.rmtree(outputs)
    manifest.pop("comparison", None)
    manifest.pop("outputs", None)


def measure(
    name: str,
    bodies: list[str] | None,
    method: str | None,
    *,
    use_saved: bool = False,
) -> None:
    """Measure a selected run, or copy the fixed saved measurements."""

    run_dir = _run_dir(name)
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = _manifest(run_dir, name)
    selected_bodies = _selected_bodies(bodies, manifest)
    selected_methods = _selected_methods(method, manifest)
    changed = False

    for body in selected_bodies:
        for selected_method in selected_methods:
            key = selected_method.lower()
            output = run_dir / "measurements" / body / f"{key}.json"
            status_output = output.with_name(f"{key}_run_status.json")
            result_record = manifest["results"].setdefault(body, {}).get(selected_method)
            if result_record and output.is_file():
                print(f"Measure: keeping existing {body} {selected_method} result")
                continue

            output.parent.mkdir(parents=True, exist_ok=True)
            if use_saved:
                defaults = APBT_PATHS[body] if selected_method == "APBT" else CIRCLE_PATHS[body]
                source = project_path(defaults.output)
                shutil.copyfile(source, output)
                source_status = source.with_name(f"{source.stem}_run_status.json")
                if source_status.is_file():
                    shutil.copyfile(source_status, status_output)
                mode = "saved_release"
                print(f"Measure: copied saved {body} {selected_method} measurements")
            else:
                module_name = (
                    "scripts.run_apbt_v2_batch"
                    if selected_method == "APBT"
                    else "scripts.run_circle_batch"
                )
                runner = import_module(module_name)
                arguments = [
                    "--body",
                    body,
                    "--output",
                    str(output),
                    "--status-output",
                    str(status_output),
                ]
                if output.exists() and status_output.exists():
                    arguments.append("--resume")
                else:
                    arguments.append("--fresh")
                return_code = runner.main(arguments)
                if return_code not in (0, 2) or not output.is_file():
                    raise WorkflowError(f"{body} {selected_method} measurement failed")
                mode = "calculated"

            record = {"path": _portable(output), "mode": mode}
            if status_output.is_file():
                record["status"] = _portable(status_output)
            manifest["results"].setdefault(body, {})[selected_method] = record
            changed = True

    selected_body_set = set(manifest.get("selection", {}).get("bodies", []))
    selected_method_set = set(manifest.get("selection", {}).get("methods", []))
    manifest["selection"] = {
        "bodies": [body for body in BODIES if body in selected_body_set | set(selected_bodies)],
        "methods": [
            value for value in METHODS if value in selected_method_set | set(selected_methods)
        ],
    }
    if changed:
        _invalidate_derived(run_dir, manifest)
    _save_manifest(run_dir, manifest)
    _remove_appledouble(run_dir)
    print(f"Measure complete: {run_dir}")


def _comparison_cases() -> list[dict]:
    return list(CASES)


def _reference_matches(requested: list[str], reference: str) -> bool:
    if "all" in requested:
        return True
    normalized = re.sub(r"[^a-z]", "", reference.lower())
    for value in requested:
        token = re.sub(r"[^a-z]", "", value.lower())
        if token == normalized or (token == "hirata" and normalized.endswith("hirata")):
            return True
    return False


def compare(
    name: str,
    bodies: list[str] | None,
    method: str | None,
    references: list[str] | None,
    metric: str | None,
) -> None:
    """Compare selected measurements against selected literature sources."""

    run_dir = _run_dir(name)
    manifest = _manifest(run_dir, name)
    if not (run_dir / "run.json").is_file():
        raise WorkflowError(f"Run {name!r} does not exist. Run the measure stage first.")
    selected_bodies = _selected_bodies(bodies, manifest)
    selected_methods = _selected_methods(method, manifest)
    requested_references = [value.lower() for value in (references or ["all"])]
    selected_metrics = {"diameter", "depth"} if metric in (None, "both") else {metric}

    raw: dict[str, dict[str, dict]] = {}
    summaries: dict[str, dict[str, dict[str, list[dict]]]] = {}
    for body in selected_bodies:
        raw[body] = {}
        summaries[body] = {}
        for selected_method in selected_methods:
            data = json_read(_result_path(manifest, body, selected_method))
            raw[body][selected_method] = data
            summaries[body][selected_method] = {"all_valid": summary_rows(body, selected_method, data)}

        if set(selected_methods) == set(METHODS):
            common = set(result_index(raw[body]["APBT"])[0]) & set(
                result_index(raw[body]["CIRCLE"])[0]
            )
            for selected_method in selected_methods:
                data = raw[body][selected_method]
                paired = filter_results(
                    data,
                    craters={Path(value).stem.lower() for value in data},
                    identities=common,
                )
                summaries[body][selected_method]["paired_common_seeds"] = summary_rows(
                    body, selected_method, paired
                )

    rows = []
    for case in _comparison_cases():
        body = case["asteroid"]
        selected_method = case["method"]
        if body not in selected_bodies or selected_method not in selected_methods:
            continue
        if case["metric"] not in selected_metrics:
            continue
        if not _reference_matches(requested_references, case["reference"]):
            continue
        for dataset, summary in summaries[body][selected_method].items():
            values = case_metrics(case, summary)
            rows.append(
                {
                    "body": body,
                    "method": selected_method,
                    "metric": case["metric"],
                    "reference": case["reference"],
                    "dataset": dataset,
                    **values,
                }
            )

    if not rows:
        raise WorkflowError(
            "No literature comparison matches that body, method, reference, and metric."
        )

    comparison_path = run_dir / "comparison.json"
    _write_json(
        comparison_path,
        {
            "selection": {
                "bodies": selected_bodies,
                "methods": selected_methods,
                "references": references or ["all"],
                "metrics": sorted(selected_metrics),
            },
            "comparisons": rows,
        },
    )
    manifest["comparison"] = {
        "path": _portable(comparison_path),
        "cases": len(rows),
    }
    manifest.pop("outputs", None)
    output_dir = run_dir / "outputs"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    _save_manifest(run_dir, manifest)
    _remove_appledouble(run_dir)

    print("Comparison complete:")
    for row in rows:
        if row["dataset"] == "all_valid":
            mae = "n/a" if row["mean_absolute_error"] is None else f"{row['mean_absolute_error']:.3f} m"
            print(
                f"- {row['body']} {row['method']} {row['metric']} vs "
                f"{row['reference']}: n={row['n']}, MAE={mae}"
            )


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def _make_figures(output_dir: Path, rows: list[dict]) -> int:
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "apbt_workflow_matplotlib")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups: dict[tuple[str, str, str], list[dict]] = {}
    for row in rows:
        if row["dataset"] == "all_valid" and row["n"]:
            key = (row["body"], row["metric"], row["reference"])
            groups.setdefault(key, []).append(row)

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    for (body, metric, reference), group in groups.items():
        fig, axis = plt.subplots(figsize=(7.2, 6.2))
        all_values = []
        for row in group:
            literature = [detail["literature"] for detail in row["details"]]
            measured = [detail["measured"] for detail in row["details"]]
            all_values.extend(literature)
            all_values.extend(measured)
            axis.scatter(literature, measured, s=42, alpha=0.75, label=row["method"])
        low = min(all_values)
        high = max(all_values)
        padding = (high - low) * 0.05 or max(abs(high) * 0.05, 1.0)
        axis.plot([low - padding, high + padding], [low - padding, high + padding], "k--")
        axis.set_xlabel(f"{reference} {metric} (m)")
        axis.set_ylabel(f"Measured {metric} (m)")
        axis.set_title(f"{body}: {metric} comparison")
        axis.grid(alpha=0.2)
        axis.legend()
        fig.tight_layout()
        fig.savefig(
            figure_dir / f"{_slug(body)}_{_slug(reference)}_{metric}.png",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(fig)
    return len(groups)


def outputs(name: str, bodies: list[str] | None, method: str | None) -> None:
    """Create final CSV tables and figures for a compared run."""

    run_dir = _run_dir(name)
    manifest = _manifest(run_dir, name)
    comparison_record = manifest.get("comparison")
    if not comparison_record:
        raise WorkflowError("No comparison found. Run the compare stage first.")
    comparison = _read_json(_stored_path(comparison_record["path"]))
    selected_bodies = _selected_bodies(
        bodies, {"selection": comparison["selection"]}
    )
    selected_methods = _selected_methods(
        method, {"selection": comparison["selection"]}
    )

    output_dir = run_dir / "outputs"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    summary_count = 0
    for body in selected_bodies:
        for selected_method in selected_methods:
            data = json_read(_result_path(manifest, body, selected_method))
            rows = summary_rows(body, selected_method, data)
            write_summary(
                rows,
                output_dir / "summaries" / body / f"{selected_method.lower()}.csv",
            )
            summary_count += 1

    comparison_rows = [
        row
        for row in comparison["comparisons"]
        if row["body"] in selected_bodies and row["method"] in selected_methods
    ]
    metric_rows = [{key: value for key, value in row.items() if key != "details"} for row in comparison_rows]
    detail_rows = []
    for row in comparison_rows:
        base = {key: value for key, value in row.items() if key != "details"}
        for detail in row["details"]:
            detail_rows.append({**base, **detail})
    csv_write(output_dir / "literature_metrics.csv", metric_rows)
    csv_write(output_dir / "literature_per_crater.csv", detail_rows)
    figure_count = _make_figures(output_dir, comparison_rows)

    report = [
        f"# Run {name}",
        "",
        f"- Bodies: {', '.join(selected_bodies)}",
        f"- Methods: {', '.join(selected_methods)}",
        f"- Summary tables: {summary_count}",
        f"- Literature cases: {len(metric_rows)}",
        f"- Figures: {figure_count}",
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(report))
    manifest["outputs"] = {
        "path": _portable(output_dir),
        "summary_tables": summary_count,
        "figures": figure_count,
    }
    _save_manifest(run_dir, manifest)
    _remove_appledouble(run_dir)
    print(f"Outputs complete: {output_dir} ({summary_count} summaries, {figure_count} figures)")


def publication_analysis(name: str, bodies: list[str] | None, method: str | None) -> None:
    """Build the legacy publication tables and figures from a complete run."""

    run_dir = _run_dir(name)
    manifest = _manifest(run_dir, name)
    selected_bodies = _selected_bodies(bodies, manifest)
    selected_methods = _selected_methods(method, manifest)
    if set(selected_bodies) != set(BODIES) or set(selected_methods) != set(METHODS):
        print("Publication analysis skipped: select all bodies and both methods.")
        return

    for body in BODIES:
        for selected_method in METHODS:
            source = run_dir / "measurements" / body / f"{selected_method.lower()}.json"
            defaults = APBT_PATHS[body] if selected_method == "APBT" else CIRCLE_PATHS[body]
            destination = project_path(defaults.output)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)

    from pipeline import analysis

    analysis.main()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Measure, compare with literature, and create selected outputs."
    )
    commands = parser.add_subparsers(dest="command", required=True)

    def selection(command, *, literature: bool = False) -> None:
        command.add_argument("--run", required=True, help="Name under runs/")
        command.add_argument("--body", nargs="+", type=_body, metavar="BODY")
        command.add_argument("--method", choices=("apbt", "circle", "both"))
        if literature:
            command.add_argument("--reference", nargs="+", default=["all"])
            command.add_argument("--metric", choices=("diameter", "depth", "both"), default="both")

    measure_parser = commands.add_parser("measure", help="Calculate selected measurements")
    selection(measure_parser)
    measure_parser.add_argument(
        "--use-saved", action="store_true", help="Copy the fixed saved results instead of recalculating"
    )

    compare_parser = commands.add_parser("compare", help="Compare a run with literature")
    selection(compare_parser, literature=True)

    output_parser = commands.add_parser("outputs", help="Create CSV and PNG outputs")
    selection(output_parser)

    all_parser = commands.add_parser("all", help="Run measure, compare, and outputs")
    selection(all_parser, literature=True)
    all_parser.add_argument(
        "--use-saved", action="store_true", help="Use fixed saved measurements for a fast run"
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "measure":
            measure(args.run, args.body, args.method, use_saved=args.use_saved)
        elif args.command == "compare":
            compare(args.run, args.body, args.method, args.reference, args.metric)
        elif args.command == "outputs":
            outputs(args.run, args.body, args.method)
        elif args.command == "all":
            measure(args.run, args.body, args.method, use_saved=args.use_saved)
            compare(args.run, args.body, args.method, args.reference, args.metric)
            outputs(args.run, args.body, args.method)
            publication_analysis(args.run, args.body, args.method)
    except (WorkflowError, FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))
    return 0
