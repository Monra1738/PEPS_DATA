"""Run the curvature-guided CIRCLE baseline for an asteroid crater set."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pyvista as pv

from methods.circle import find_crater_rim_from_seeds, project_circle_with_rays
from scripts.batch_config import CIRCLE_PATHS, portable_path, project_path
from utils.mesh_manipulation import (
    curv_copy_mesh,
    force_meters,
    frac_neighbors_from_radius,
    return_depth,
    rotate_mesh,
    scale_from_source_unit_to_meters,
    signal_choice,
)


class InvalidCircleMeasurementError(ValueError):
    """Raised when derived CIRCLE rim/depth geometry is not physically valid."""


def validate_circle_measurements(
    *,
    depth_med: float,
    depth: float,
    diameter_max: float,
) -> None:
    """Reject non-finite or non-positive CIRCLE geometry before recording it."""

    measurements = {
        "Depth_med": float(depth_med),
        "Depth": float(depth),
        "Diameter_max": float(diameter_max),
    }
    nonfinite = {
        name: value for name, value in measurements.items() if not np.isfinite(value)
    }
    if nonfinite:
        raise InvalidCircleMeasurementError(
            f"Non-finite CIRCLE rim/depth measurement: {nonfinite}"
        )
    nonpositive = {
        name: value for name, value in measurements.items() if value <= 0.0
    }
    if nonpositive:
        raise InvalidCircleMeasurementError(
            "Non-positive CIRCLE measurement indicates invalid rim/floor geometry: "
            f"{nonpositive}"
        )


def _normalize_key(path_like: str) -> str:
    return Path(path_like).stem.lower()


def _normalize_entries(value):
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    if isinstance(value, dict):
        return [value]
    return []


def _build_lookup(config_dict: dict):
    return {
        _normalize_key(full_path): _normalize_entries(entry)
        for full_path, entry in config_dict.items()
    }


def _pick_entry(entries, idx: int):
    if not entries:
        return None
    if idx < len(entries):
        return entries[idx]
    return entries[-1]


def load_ev_param(crater_name: str, seed_idx: int, cfg_lookup: dict):
    selected_cfg = _pick_entry(cfg_lookup.get(_normalize_key(crater_name), []), seed_idx)
    if selected_cfg is None:
        raise ValueError(f"No configuration found for crater {crater_name!r}")

    if "TAUBIN_N_ITER" in selected_cfg:
        return (
            int(selected_cfg.get("TAUBIN_N_ITER", 1)),
            float(selected_cfg.get("TAUBIN_PASS_BAND", 0.1)),
            int(selected_cfg.get("POOL_RINGS", 1)),
            str(selected_cfg.get("POOL_METHOD", "median")),
            str(selected_cfg.get("SIGNAL_CHOICE", "abs_k1")),
            bool(selected_cfg.get("INVERSE", False)),
            selected_cfg.get("SCALE_HINT", "m"),
        )

    smoothing_cfg = selected_cfg["Smoothing: "]
    pooling_cfg = selected_cfg["Pooling : "]
    return (
        int(smoothing_cfg["n_iter"]),
        float(smoothing_cfg.get("pass_band", 0.1)),
        int(pooling_cfg["number_rings"]),
        str(pooling_cfg.get("pool_method", "median")),
        str(selected_cfg["Signal"]),
        bool(selected_cfg.get("Mesh_inverse", False)),
        selected_cfg.get("SCALE_HINT", "m"),
    )


def _write_json_atomic(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as fh:
        json.dump(value, fh, indent=2, default=float)
        fh.write("\n")
    os.replace(temporary, path)


def _load_output(path: Path, *, fresh: bool, overwrite: bool) -> dict:
    if fresh:
        if path.exists() and not overwrite:
            raise FileExistsError(
                f"Fresh output already exists: {path}. Use --overwrite only for "
                "an explicitly disposable verification output."
            )
        return {}
    if not path.exists():
        return {}
    try:
        with path.open() as fh:
            return json.load(fh)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Existing output is not valid JSON: {path}") from exc


def _seed_identity(crater: str, seed_points) -> tuple[str, tuple[int, ...]]:
    return _normalize_key(crater), tuple(int(value) for value in seed_points)


def _prepare_resume_attempts(
    db: dict,
    attempts: list[dict],
) -> tuple[list[dict], set[tuple[str, tuple[int, ...]]]]:
    """Validate saved progress and return attempts plus identities to skip."""

    attempt_index = {}
    retained_attempts = []
    for attempt in attempts:
        identity = _seed_identity(attempt["crater"], attempt["seed_points"])
        if identity in attempt_index:
            raise ValueError(f"Duplicate saved CIRCLE attempt identity: {identity}")
        if attempt.get("status") == "interrupted":
            continue
        attempt_index[identity] = attempt
        retained_attempts.append(attempt)

    result_index = {}
    for crater, results in db.items():
        for result in results:
            identity = _seed_identity(crater, result["seed_points"])
            if identity in result_index:
                raise ValueError(f"Duplicate saved CIRCLE result identity: {identity}")
            result_index[identity] = result

    for identity, result in result_index.items():
        attempt = attempt_index.get(identity)
        if attempt is not None and attempt.get("status") != "success":
            raise ValueError(
                "Saved CIRCLE result conflicts with a non-success status for "
                f"{identity}"
            )
        if attempt is None:
            crater, seed_points = identity
            attempt = {
                "crater": crater,
                "seed_points": list(seed_points),
                "status": "success",
                "diameter_max": result.get("Diameter_max"),
                "depth_med": result.get("Depth_med"),
                "depth": result.get("Depth"),
                "recovered_from_result": True,
            }
            attempt_index[identity] = attempt
            retained_attempts.append(attempt)

    for identity, attempt in attempt_index.items():
        if attempt.get("status") == "success" and identity not in result_index:
            raise ValueError(
                f"Saved successful CIRCLE status has no matching result: {identity}"
            )

    return retained_attempts, set(attempt_index)


def _scale_mesh(mesh: pv.PolyData, source_unit: str | None) -> pv.PolyData:
    if source_unit:
        return scale_from_source_unit_to_meters(mesh, source_unit)
    # This reproduces the final lightweight runner. It intentionally does not
    # pass the historical SCALE_HINT because those labels used inconsistent
    # multipliers; all current local asteroid meshes are stored in kilometres.
    return force_meters(mesh)


def run_batch(
    *,
    mesh_dir: Path,
    seeds_path: Path,
    config_path: Path,
    output_path: Path,
    status_path: Path,
    fresh: bool,
    overwrite: bool,
    resume: bool = False,
    source_unit: str | None,
    crater_filters: set[str] | None = None,
) -> dict:
    if fresh and resume:
        raise ValueError("fresh and resume modes are mutually exclusive")
    with seeds_path.open() as fh:
        seed_data = json.load(fh)
    with config_path.open() as fh:
        config_data = json.load(fh)

    cfg_lookup = _build_lookup(config_data)
    vtk_index = {path.stem.lower(): path for path in mesh_dir.glob("*.vtk")}
    db = _load_output(output_path, fresh=fresh, overwrite=overwrite)

    status_metadata = {
        "method": "CIRCLE",
        "mesh_dir": portable_path(mesh_dir),
        "seeds": portable_path(seeds_path),
        "config": portable_path(config_path),
        "output": portable_path(output_path),
        "source_unit_mode": source_unit or "legacy force_meters",
        "started_unix": time.time(),
        "attempts": [],
    }
    completed_identities = set()
    if resume and status_path.exists():
        with status_path.open() as fh:
            saved_status = json.load(fh)
        comparable_fields = (
            "method",
            "mesh_dir",
            "seeds",
            "config",
            "output",
            "source_unit_mode",
        )
        mismatches = {
            field: {"saved": saved_status.get(field), "requested": status_metadata.get(field)}
            for field in comparable_fields
            if saved_status.get(field) != status_metadata.get(field)
        }
        if mismatches:
            raise ValueError(
                f"Resume metadata does not match saved CIRCLE run: {mismatches}"
            )
        attempts, completed_identities = _prepare_resume_attempts(
            db,
            saved_status.get("attempts", []),
        )
        status = {
            **status_metadata,
            "started_unix": saved_status.get("started_unix", time.time()),
            "attempts": attempts,
            "resume_events_unix": [
                *saved_status.get("resume_events_unix", []),
                time.time(),
            ],
        }
    elif resume:
        attempts, completed_identities = _prepare_resume_attempts(db, [])
        status = {
            **status_metadata,
            "attempts": attempts,
            "resume_events_unix": [time.time()],
        }
    else:
        status = status_metadata
    _write_json_atomic(status_path, status)

    for seed_key, seed_groups in seed_data.items():
        crater_name = Path(seed_key).stem
        normalized_crater = _normalize_key(crater_name)
        if crater_filters and normalized_crater not in crater_filters:
            continue
        mesh_path = vtk_index.get(normalized_crater)
        if mesh_path is None:
            for seed_idx, seed in enumerate(seed_groups):
                status["attempts"].append({
                    "crater": crater_name,
                    "seed_index": seed_idx,
                    "seed_points": seed,
                    "status": "failed",
                    "reason": "mesh_not_found",
                })
            _write_json_atomic(status_path, status)
            continue

        crater_results = db.setdefault(crater_name, [])
        for seed_idx, seed in enumerate(seed_groups):
            identity = _seed_identity(crater_name, seed)
            if identity in completed_identities:
                print(
                    f"CIRCLE resume: skipping saved {crater_name}#{seed_idx}",
                    flush=True,
                )
                continue
            started = time.perf_counter()
            attempt = {
                "crater": crater_name,
                "seed_index": seed_idx,
                "seed_points": seed,
                "mesh": portable_path(mesh_path),
            }
            try:
                (
                    taubin_iterations,
                    taubin_pass_band,
                    pool_rings,
                    pool_method,
                    signal_name,
                    inverse,
                    configured_scale_hint,
                ) = load_ev_param(crater_name, seed_idx, cfg_lookup)

                points_circle = 200
                ray_length = 1e3
                curvature_low = 0.70
                curvature_high = 0.99
                rim_band = 0.25

                mesh = rotate_mesh(pv.read(mesh_path), inverse)
                mesh = _scale_mesh(mesh, source_unit)
                _, _, mean_curvature, _, k1, k2 = curv_copy_mesh(
                    mesh,
                    taubin_iterations,
                    taubin_pass_band,
                    pool_rings,
                    pool_method,
                )
                curvature_signal = signal_choice(
                    mean_curvature,
                    k1,
                    k2,
                    s_choice=signal_name,
                )
                center_xy, radius, _, diameter, rim_indices = find_crater_rim_from_seeds(
                    mesh,
                    seed,
                    curvature_signal,
                    rim_band,
                    curvature_low,
                    curvature_high,
                )
                if center_xy is None:
                    raise ValueError("insufficient_curvature_candidates")

                rim_points = np.asarray(mesh.points[rim_indices])
                if rim_points.ndim != 2 or rim_points.shape[0] == 0:
                    raise ValueError("empty_curvature_rim")
                rim_z_reference = np.median(rim_points[:, 2])
                theta = np.linspace(0, 2 * np.pi, points_circle, endpoint=True)
                circle_xy = np.column_stack([
                    center_xy[0] + radius * np.cos(theta),
                    center_xy[1] + radius * np.sin(theta),
                ])
                circle_xyz = np.column_stack([
                    circle_xy,
                    np.full_like(theta, rim_z_reference),
                ])
                final_rim = np.asarray(
                    project_circle_with_rays(
                        mesh,
                        circle_xyz,
                        ray_length,
                        direction=(0, 0, -1),
                    )
                )
                if final_rim.ndim != 2 or final_rim.shape[0] == 0:
                    raise ValueError("empty_projected_rim")

                fraction = frac_neighbors_from_radius(radius)
                depth_result = return_depth(mesh, final_rim, center_xy, fraction)
                if depth_result[0] is None:
                    raise ValueError("depth_calculation_failed")
                (
                    depth_med,
                    depth,
                    depth_p40,
                    depth_p60,
                    rim_level,
                    floor_level,
                    err_minus,
                    err_plus,
                    fraction,
                ) = depth_result
                diameter_max = radius * 2
                validate_circle_measurements(
                    depth_med=depth_med,
                    depth=depth,
                    diameter_max=diameter_max,
                )

                metadata = {
                    "seed_points": np.asarray(seed).tolist(),
                    "inverse": inverse,
                    "POINTS_CIRCLE": points_circle,
                    "RAY_LENGTH": ray_length,
                    "Diameter_max": diameter_max,
                    "Depth_med": depth_med,
                    "Depth": depth,
                    "Depth_p40": float(depth_p40),
                    "Depth_p60": float(depth_p60),
                    "Rim_level": rim_level,
                    "Floor_level": floor_level,
                    "Err_minus": err_minus,
                    "Err_plus": err_plus,
                    "Frac_neighbors": fraction,
                }
                crater_results.append(metadata)
                completed_identities.add(identity)
                _write_json_atomic(output_path, db)
                attempt.update({
                    "status": "success",
                    "diameter_max": float(diameter),
                    "depth_med": float(depth_med),
                    "configured_scale_hint": configured_scale_hint,
                })
            except KeyboardInterrupt:
                attempt.update({
                    "status": "interrupted",
                    "reason": "KeyboardInterrupt",
                    "message": "Run interrupted by operator",
                })
                raise
            except Exception as exc:
                message = str(exc)
                reason = (
                    type(exc).__name__
                    if isinstance(exc, InvalidCircleMeasurementError)
                    else (message.split(":", 1)[0] if message else type(exc).__name__)
                )
                attempt.update({
                    "status": "failed",
                    "reason": reason,
                    "exception_type": type(exc).__name__,
                    "message": message,
                })
                print(
                    f"ERROR CIRCLE {crater_name} seed #{seed_idx} {seed}: {exc}",
                    flush=True,
                )
                traceback.print_exc()
            finally:
                attempt["elapsed_seconds"] = time.perf_counter() - started
                status["attempts"].append(attempt)
                success_count = sum(
                    row.get("status") == "success" for row in status["attempts"]
                )
                print(
                    f"CIRCLE progress: {len(status['attempts'])} attempted, "
                    f"{success_count} successful; latest={crater_name}#{seed_idx}",
                    flush=True,
                )
                _write_json_atomic(status_path, status)

    status["finished_unix"] = time.time()
    status["attempted"] = len(status["attempts"])
    status["successful"] = sum(
        row.get("status") == "success" for row in status["attempts"]
    )
    status["failed"] = status["attempted"] - status["successful"]
    _write_json_atomic(status_path, status)
    if not output_path.exists():
        _write_json_atomic(output_path, db)
    return status


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--body", choices=sorted(CIRCLE_PATHS), default="Bennu")
    parser.add_argument("--mesh-dir", type=Path)
    parser.add_argument("--seeds", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--status-output", type=Path)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume an interrupted output by skipping validated saved seed identities.",
    )
    parser.add_argument(
        "--source-unit",
        choices=("m", "km", "cm", "mm"),
        help="Use explicit physical units instead of historical scaling.",
    )
    parser.add_argument(
        "--crater",
        action="append",
        default=[],
        help="Limit the run to a crater name; may be supplied more than once.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    defaults = CIRCLE_PATHS[args.body]
    mesh_dir = project_path(args.mesh_dir or defaults.mesh_dir)
    seeds = project_path(args.seeds or defaults.seeds)
    config = project_path(args.config or defaults.config)
    output = project_path(args.output or defaults.output)
    status_output = project_path(
        args.status_output or output.with_name(output.stem + "_run_status.json")
    )
    crater_filters = {_normalize_key(value) for value in args.crater} or None

    status = run_batch(
        mesh_dir=mesh_dir,
        seeds_path=seeds,
        config_path=config,
        output_path=output,
        status_path=status_output,
        fresh=args.fresh,
        overwrite=args.overwrite,
        resume=args.resume,
        source_unit=args.source_unit,
        crater_filters=crater_filters,
    )
    print(
        f"CIRCLE complete: {status['successful']}/{status['attempted']} successful; "
        f"results={output}; status={status_output}"
    )
    return 0 if status["failed"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
