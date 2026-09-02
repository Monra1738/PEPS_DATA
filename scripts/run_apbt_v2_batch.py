"""Run APBT v2 for a configured asteroid crater set.

The original research script executed immediately at import time and was
hard-coded to Bennu.  This entry point retains the numerical workflow while
making inputs and outputs explicit, allowing fresh verification runs without
touching the canonical result JSONs.
"""

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

from methods.apbt_v2 import (
    _based_diffusion_matrix,
    _find_center_radius_from_segement_rim,
    apply_gradient_diffusion,
    calculate_rim_diameter,
    detrain_mesh,
    estimate_the_best_center,
    find_initial_center,
    find_inside_pts,
    find_the_best_seed_points,
    frac_neighbors_from_radius,
    get_best_rim_pts_per_bin,
    grow_segments_unitl_anchor,
    project_polyline_with_rays,
    return_poly_rim,
    smooth_rim_angular,
    upsample_rim_by_angle,
)
from scripts.batch_config import APBT_PATHS, portable_path, project_path
from utils.mesh_manipulation import (
    adj_vertex_list,
    compute_reference_normals,
    force_meters,
    return_depth,
    rotate_mesh_and_normals,
    scale_from_source_unit_to_meters,
)


HISTORICAL_CENTER_N_ITER = 60
HISTORICAL_CENTER_WEIGHT = 0.5
CENTER_PARAM_MODES = ("historical", "configured")


class InvalidAPBTMeasurementError(ValueError):
    """Raised when the derived rim/depth geometry is not physically valid."""


def validate_apbt_measurements(
    *,
    depth_med: float,
    depth: float,
    diameter_ew: float,
    diameter_ns: float,
    diameter_max: float,
) -> None:
    """Reject non-finite or non-positive APBT geometry before it is recorded."""

    measurements = {
        "Depth_med": float(depth_med),
        "Depth": float(depth),
        "Diameter_ew": float(diameter_ew),
        "Diameter_ns": float(diameter_ns),
        "Diameter_max": float(diameter_max),
    }
    nonfinite = {
        name: value for name, value in measurements.items() if not np.isfinite(value)
    }
    if nonfinite:
        raise InvalidAPBTMeasurementError(
            f"Non-finite APBT rim/depth measurement: {nonfinite}"
        )
    nonpositive_depths = {
        name: measurements[name]
        for name in ("Depth_med", "Depth")
        if measurements[name] <= 0.0
    }
    if nonpositive_depths:
        raise InvalidAPBTMeasurementError(
            "Non-positive depth indicates invalid APBT rim/floor geometry: "
            f"{nonpositive_depths}"
        )
    nonpositive_diameters = {
        name: measurements[name]
        for name in ("Diameter_ew", "Diameter_ns", "Diameter_max")
        if measurements[name] <= 0.0
    }
    if nonpositive_diameters:
        raise InvalidAPBTMeasurementError(
            "Non-positive diameter indicates invalid APBT rim geometry: "
            f"{nonpositive_diameters}"
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
    lookup = {}
    for full_path, entry in config_dict.items():
        lookup[_normalize_key(full_path)] = _normalize_entries(entry)
    return lookup


def _pick_entry(entries, _index: int):
    """Retain the historical APBT behavior: use the final crater config."""

    if not entries:
        return None
    return entries[-1]


def load_ev_param(
    crater_name: str,
    seed_idx: int,
    cfg_lookup: dict,
    inverse_lookup: dict | None = None,
):
    key = _normalize_key(crater_name)
    selected_cfg = _pick_entry(cfg_lookup.get(key, []), seed_idx)

    if not selected_cfg:
        raise ValueError(f"No configuration found for crater {crater_name!r}")

    n_iter_1 = int(selected_cfg["n_iter_1"])
    n_iter_2 = int(selected_cfg["n_iter_2"])
    n_iter_3 = int(selected_cfg["n_iter_3"])
    top_percent_1 = float(selected_cfg["top_percent_1"])
    top_percent_2 = float(selected_cfg["top_percent_2"])
    things_weight = float(selected_cfg["things_weight"])
    use_detrain = bool(selected_cfg["USE_DETRAIN"])

    inverse = selected_cfg.get("INVERSE")
    if inverse is None and inverse_lookup is not None:
        inv_cfg = _pick_entry(inverse_lookup.get(key, []), seed_idx)
        if inv_cfg is not None:
            inverse = inv_cfg.get("INVERSE")
    if inverse is None:
        raise ValueError(f"No INVERSE flag found for crater {crater_name!r}")

    return (
        n_iter_1,
        n_iter_2,
        n_iter_3,
        top_percent_1,
        top_percent_2,
        things_weight,
        use_detrain,
        bool(inverse),
    )


def resolve_center_parameters(
    configured_n_iter: int,
    configured_weight: float,
    *,
    n_iter_mode: str,
    weight_mode: str,
) -> tuple[int, float]:
    """Resolve independently controlled center-refinement parameters."""

    if n_iter_mode not in CENTER_PARAM_MODES:
        raise ValueError(
            f"Unsupported center iteration mode {n_iter_mode!r}; "
            f"expected one of {CENTER_PARAM_MODES}"
        )
    if weight_mode not in CENTER_PARAM_MODES:
        raise ValueError(
            f"Unsupported center weight mode {weight_mode!r}; "
            f"expected one of {CENTER_PARAM_MODES}"
        )
    effective_n_iter = (
        int(configured_n_iter)
        if n_iter_mode == "configured"
        else HISTORICAL_CENTER_N_ITER
    )
    effective_weight = (
        float(configured_weight)
        if weight_mode == "configured"
        else HISTORICAL_CENTER_WEIGHT
    )
    return effective_n_iter, effective_weight


def resolve_center_parameter_modes(
    *,
    use_configured_center_params: bool,
    center_n_iter_mode: str | None,
    center_weight_mode: str | None,
) -> tuple[str, str]:
    """Resolve the legacy combined flag and the independent mode flags."""

    if use_configured_center_params:
        if center_n_iter_mode not in (None, "configured"):
            raise ValueError(
                "--use-configured-center-params conflicts with an explicit "
                "historical center iteration mode"
            )
        if center_weight_mode not in (None, "configured"):
            raise ValueError(
                "--use-configured-center-params conflicts with an explicit "
                "historical center weight mode"
            )
        center_n_iter_mode = "configured"
        center_weight_mode = "configured"
    return center_n_iter_mode or "historical", center_weight_mode or "historical"


def _sanitize_seed_points(seed_points, n_points: int):
    valid = []
    invalid = []
    for value in seed_points:
        try:
            idx = int(value)
        except (TypeError, ValueError):
            invalid.append(value)
            continue
        if 0 <= idx < n_points:
            valid.append(idx)
        else:
            invalid.append(idx)
    return valid, invalid


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
            raise ValueError(f"Duplicate saved APBT attempt identity: {identity}")
        if attempt.get("status") == "interrupted":
            continue
        attempt_index[identity] = attempt
        retained_attempts.append(attempt)

    result_index = {}
    for crater, results in db.items():
        for result in results:
            identity = _seed_identity(crater, result["seed_points"])
            if identity in result_index:
                raise ValueError(f"Duplicate saved APBT result identity: {identity}")
            result_index[identity] = result

    for identity, result in result_index.items():
        attempt = attempt_index.get(identity)
        if attempt is not None and attempt.get("status") != "success":
            raise ValueError(
                "Saved APBT result conflicts with a non-success status for "
                f"{identity}"
            )
        if attempt is None:
            crater, seed_points = identity
            attempt = {
                "crater": crater,
                "seed_points": list(seed_points),
                "variant_name": result.get("PARAMETER_VARIANT"),
                "configured_n_iter_1": result.get(
                    "CONFIGURED_N_ITER_1", result.get("N_ITER_1")
                ),
                "configured_things_weight": result.get(
                    "CONFIGURED_THINGES_WEIGHT", result.get("THINGES_WEIGHT")
                ),
                "effective_n_iter_1": result.get("EFFECTIVE_N_ITER_1"),
                "effective_things_weight": result.get("EFFECTIVE_THINGES_WEIGHT"),
                "center_n_iter_mode": result.get("CENTER_N_ITER_MODE"),
                "center_weight_mode": result.get("CENTER_WEIGHT_MODE"),
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
                f"Saved successful APBT status has no matching result: {identity}"
            )

    return retained_attempts, set(attempt_index)


def _scale_mesh(mesh: pv.PolyData, source_unit: str | None) -> pv.PolyData:
    if source_unit:
        return scale_from_source_unit_to_meters(mesh, source_unit)
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
    use_configured_center_params: bool = False,
    center_n_iter_mode: str | None = None,
    center_weight_mode: str | None = None,
    variant_name: str | None = None,
    crater_filters: set[str] | None = None,
) -> dict:
    if fresh and resume:
        raise ValueError("fresh and resume modes are mutually exclusive")
    center_n_iter_mode, center_weight_mode = resolve_center_parameter_modes(
        use_configured_center_params=use_configured_center_params,
        center_n_iter_mode=center_n_iter_mode,
        center_weight_mode=center_weight_mode,
    )
    # Validate the modes before creating any output file.
    resolve_center_parameters(
        HISTORICAL_CENTER_N_ITER,
        HISTORICAL_CENTER_WEIGHT,
        n_iter_mode=center_n_iter_mode,
        weight_mode=center_weight_mode,
    )

    with seeds_path.open() as fh:
        seed_data = json.load(fh)
    with config_path.open() as fh:
        config_data = json.load(fh)

    cfg_lookup = _build_lookup(config_data)
    inverse_lookup = cfg_lookup
    vtk_index = {path.stem.lower(): path for path in mesh_dir.glob("*.vtk")}
    db = _load_output(output_path, fresh=fresh, overwrite=overwrite)

    status_metadata = {
        "method": "APBT v2",
        "mesh_dir": portable_path(mesh_dir),
        "seeds": portable_path(seeds_path),
        "config": portable_path(config_path),
        "output": portable_path(output_path),
        "source_unit_mode": source_unit or "legacy force_meters",
        "configured_center_params": (
            center_n_iter_mode == "configured"
            and center_weight_mode == "configured"
        ),
        "center_n_iter_mode": center_n_iter_mode,
        "center_weight_mode": center_weight_mode,
        "variant_name": variant_name,
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
            "center_n_iter_mode",
            "center_weight_mode",
            "variant_name",
        )
        mismatches = {
            field: {"saved": saved_status.get(field), "requested": status_metadata.get(field)}
            for field in comparable_fields
            if saved_status.get(field) != status_metadata.get(field)
        }
        if mismatches:
            raise ValueError(f"Resume metadata does not match saved APBT run: {mismatches}")
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
                    f"APBT resume: skipping saved {crater_name}#{seed_idx}",
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
                    n_iter_1,
                    n_iter_2,
                    n_iter_3,
                    top_percent_1,
                    top_percent_2,
                    things_weight,
                    use_detrain,
                    inverse,
                ) = load_ev_param(crater_name, seed_idx, cfg_lookup, inverse_lookup)
                effective_n_iter_1, effective_things_weight = resolve_center_parameters(
                    n_iter_1,
                    things_weight,
                    n_iter_mode=center_n_iter_mode,
                    weight_mode=center_weight_mode,
                )
                attempt.update({
                    "variant_name": variant_name,
                    "configured_n_iter_1": n_iter_1,
                    "configured_things_weight": things_weight,
                    "effective_n_iter_1": effective_n_iter_1,
                    "effective_things_weight": effective_things_weight,
                    "center_n_iter_mode": center_n_iter_mode,
                    "center_weight_mode": center_weight_mode,
                })

                points_circle = 180
                ray_length = 200

                mesh = _scale_mesh(pv.read(mesh_path), source_unit)
                valid_seeds, invalid_seeds = _sanitize_seed_points(seed, mesh.n_points)
                if invalid_seeds:
                    raise ValueError(
                        f"Invalid seed indices {invalid_seeds}; mesh has {mesh.n_points} points"
                    )
                if len(valid_seeds) < 3:
                    raise ValueError(
                        f"Need at least 3 valid seed indices, got {len(valid_seeds)}"
                    )

                normals = compute_reference_normals(mesh)
                mesh, normals = rotate_mesh_and_normals(mesh, normals, inverse=inverse)
                adjacency = adj_vertex_list(mesh)

                working_mesh = detrain_mesh(mesh) if use_detrain else mesh
                matrix = _based_diffusion_matrix(working_mesh, adjacency)
                best_center = estimate_the_best_center(
                    working_mesh,
                    normals,
                    matrix,
                    top_percent_1,
                    N_ITER_1=effective_n_iter_1,
                    THINGES_WEIGHT=effective_things_weight,
                )
                working_mesh = return_poly_rim(
                    working_mesh,
                    normals,
                    matrix,
                    best_center,
                    smooth_iter=n_iter_2,
                    top_percent_2=top_percent_2,
                )
                working_mesh, score_field = apply_gradient_diffusion(
                    working_mesh,
                    matrix,
                    n_iter=n_iter_3,
                )

                projection_mesh = working_mesh.copy()
                analysis_mesh = working_mesh.copy()
                best_seeds = find_the_best_seed_points(
                    analysis_mesh,
                    np.asarray(valid_seeds, dtype=int),
                    score_field,
                    adjacency,
                    n_rings=0,
                )
                center, radius, seed_shape = find_initial_center(analysis_mesh, best_seeds)
                low_points = find_inside_pts(analysis_mesh, center, radius, seed_shape)
                segment_rim = grow_segments_unitl_anchor(
                    analysis_mesh,
                    low_points,
                    best_seeds,
                    adjacency,
                )
                if np.asarray(segment_rim).size == 0:
                    raise ValueError("Topology growth returned an empty rim")

                center, radius, rim_shape = _find_center_radius_from_segement_rim(
                    analysis_mesh,
                    segment_rim,
                )
                thinned_rim = get_best_rim_pts_per_bin(
                    analysis_mesh,
                    segment_rim,
                    center,
                    score_field,
                    shape_info=rim_shape,
                )
                smooth_points = smooth_rim_angular(
                    analysis_mesh,
                    thinned_rim,
                    center,
                    score_field,
                    shape_info=rim_shape,
                )
                smooth_center = center.copy()
                smooth_center[2] = smooth_points[:, 2].mean()
                dense_points = upsample_rim_by_angle(
                    smooth_points,
                    smooth_center,
                    n_points=points_circle,
                    shape_info=rim_shape,
                )
                final_rim = project_polyline_with_rays(
                    projection_mesh,
                    segment_rim,
                    dense_points,
                    radius,
                    ray_length=ray_length,
                )
                if final_rim.ndim != 2 or final_rim.shape[0] == 0:
                    raise ValueError("Projected APBT rim is empty")

                fraction = frac_neighbors_from_radius(radius)
                depth_result = return_depth(
                    analysis_mesh,
                    final_rim,
                    smooth_center[:2],
                    fraction,
                )
                if depth_result[0] is None:
                    raise ValueError("Depth calculation returned no result")
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
                diameter_ew, diameter_ns, diameter_max = calculate_rim_diameter(final_rim)
                validate_apbt_measurements(
                    depth_med=depth_med,
                    depth=depth,
                    diameter_ew=diameter_ew,
                    diameter_ns=diameter_ns,
                    diameter_max=diameter_max,
                )
                metadata = {
                    "seed_points": np.asarray(seed).tolist(),
                    "FACTOR_RADIUS": 2,
                    "inverse": inverse,
                    "STD_DEV_FACTOR_ADAP": 1.0,
                    "ANGLE_SIGMA": 20.0,
                    "SMOOTH_N_ITER": 3,
                    "POINTS_CIRCLE": points_circle,
                    "RAY_LENGTH": ray_length,
                    "N_ITER_1": n_iter_1,
                    "N_ITER_2": n_iter_2,
                    "N_ITER_3": n_iter_3,
                    "TOP_PERCENT_1": top_percent_1,
                    "TOP_PERCENT_2": top_percent_2,
                    "THINGES_WEIGHT": things_weight,
                    "CONFIGURED_N_ITER_1": n_iter_1,
                    "CONFIGURED_THINGES_WEIGHT": things_weight,
                    "EFFECTIVE_N_ITER_1": effective_n_iter_1,
                    "EFFECTIVE_THINGES_WEIGHT": effective_things_weight,
                    "CENTER_N_ITER_MODE": center_n_iter_mode,
                    "CENTER_WEIGHT_MODE": center_weight_mode,
                    "PARAMETER_VARIANT": variant_name,
                    "USE_DETRAIN": use_detrain,
                    "Diameter_ew": diameter_ew,
                    "Diameter_ns": diameter_ns,
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
                    "diameter_max": float(diameter_max),
                    "depth_med": float(depth_med),
                    "depth": float(depth),
                })
            except KeyboardInterrupt:
                attempt.update({
                    "status": "interrupted",
                    "reason": "KeyboardInterrupt",
                    "message": "Run interrupted by operator",
                })
                raise
            except Exception as exc:
                attempt.update({
                    "status": "failed",
                    "reason": type(exc).__name__,
                    "message": str(exc),
                })
                print(
                    f"ERROR APBT {crater_name} seed #{seed_idx} {seed}: {exc}",
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
                    f"APBT progress: {len(status['attempts'])} attempted, "
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
    parser.add_argument("--body", choices=sorted(APBT_PATHS), default="Bennu")
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
        help="Use explicit physical units instead of the historical scaling helper.",
    )
    parser.add_argument(
        "--use-configured-center-params",
        action="store_true",
        help=(
            "Backward-compatible alias for --center-n-iter-mode configured "
            "--center-weight-mode configured."
        ),
    )
    parser.add_argument(
        "--center-n-iter-mode",
        choices=CENTER_PARAM_MODES,
        help="Use historical N_ITER_1=60 or the crater-configured value.",
    )
    parser.add_argument(
        "--center-weight-mode",
        choices=CENTER_PARAM_MODES,
        help="Use historical THINGES_WEIGHT=0.5 or the crater-configured value.",
    )
    parser.add_argument(
        "--variant-name",
        help="Optional diagnostic label stored in result and status metadata.",
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
    defaults = APBT_PATHS[args.body]
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
        use_configured_center_params=args.use_configured_center_params,
        center_n_iter_mode=args.center_n_iter_mode,
        center_weight_mode=args.center_weight_mode,
        variant_name=args.variant_name,
        crater_filters=crater_filters,
    )
    print(
        f"APBT complete: {status['successful']}/{status['attempted']} successful; "
        f"results={output}; status={status_output}"
    )
    return 0 if status["failed"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
