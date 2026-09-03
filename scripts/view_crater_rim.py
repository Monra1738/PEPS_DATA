"""Visualize one APBT crater rim.

Example:
    python scripts/view_crater_rim.py --body Bennu --crater 21 --seed 0
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

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
    get_best_rim_pts_per_bin,
    grow_segments_unitl_anchor,
    project_polyline_with_rays,
    return_poly_rim,
    smooth_rim_angular,
    upsample_rim_by_angle,
)
from scripts.batch_config import APBT_PATHS, project_path
from scripts.run_apbt_v2_batch import (
    _build_lookup,
    _scale_mesh,
    load_ev_param,
    resolve_center_parameters,
)
from utils.mesh_manipulation import (
    adj_vertex_list,
    compute_reference_normals,
    rotate_mesh_and_normals,
)


BODIES = tuple(APBT_PATHS)


def _read(path: Path):
    with path.open() as handle:
        return json.load(handle)


def _result_path(body: str, result_path: Path | None) -> Path:
    if result_path is not None:
        return project_path(result_path)
    run_candidates = sorted((project_path("runs")).glob(f"*/measurements/{body}/apbt.json"))
    if run_candidates:
        return run_candidates[-1]
    return project_path(APBT_PATHS[body].output)


def _entry(body: str, crater: str, seed_index: int, result_path: Path | None) -> dict:
    results = _read(_result_path(body, result_path))
    crater_text = str(crater)
    crater_key = crater_text if crater_text.startswith("crater_") else f"crater_{int(crater_text)}"
    entries = results.get(crater_key)
    if not entries or seed_index >= len(entries):
        raise ValueError(f"No saved APBT result for {crater_key} seed {seed_index}")
    return entries[seed_index]


def _line(points: np.ndarray) -> pv.PolyData:
    points = np.asarray(points)
    if len(points) < 2:
        raise ValueError("APBT produced fewer than two rim points")
    line = pv.PolyData(points)
    indices = np.arange(len(points), dtype=np.int64)
    line.lines = np.concatenate(([len(indices)], indices))
    return line


def _debug_value(label: str, value) -> None:
    """Print a compact description of a value returned by an APBT step."""
    if isinstance(value, np.ndarray):
        finite = value[np.isfinite(value)] if np.issubdtype(value.dtype, np.number) else []
        details = f"shape={value.shape}, dtype={value.dtype}"
        if len(finite):
            details += f", min={finite.min():.6g}, max={finite.max():.6g}"
        print(f"[visualize] {label}: ndarray({details})")
    elif hasattr(value, "n_points") and hasattr(value, "n_cells"):
        print(f"[visualize] {label}: PolyData(points={value.n_points}, cells={value.n_cells}, fields={list(value.point_data)})")
    elif hasattr(value, "shape") and hasattr(value, "nnz"):
        print(f"[visualize] {label}: sparse(shape={value.shape}, nnz={value.nnz})")
    elif isinstance(value, dict):
        print(f"[visualize] {label}: dict(keys={list(value)})")
    elif isinstance(value, (tuple, list)):
        print(f"[visualize] {label}: {type(value).__name__}(len={len(value)})")
        for index, item in enumerate(value):
            _debug_value(f"{label}[{index}]", item)
    else:
        print(f"[visualize] {label}: {value!r}")


def calculate_rim(body: str, crater: str, seed_index: int, result_path: Path | None, *, debug: bool = True):
    def report(label: str, value) -> None:
        if debug:
            _debug_value(label, value)

    defaults = APBT_PATHS[body]
    crater_name = str(crater)
    if not crater_name.startswith("crater_"):
        crater_name = f"crater_{int(crater_name)}"
    mesh_path = project_path(defaults.mesh_dir) / f"{crater_name}.vtk"
    if not mesh_path.is_file():
        raise FileNotFoundError(mesh_path)

    result = _entry(body, crater_name, seed_index, result_path)
    seed = np.asarray(result["seed_points"], dtype=int)
    mesh = _scale_mesh(pv.read(mesh_path), None)
    report("loaded mesh", mesh)
    report("saved APBT result", result)
    report("seed points", seed)
    if np.any(seed < 0) or np.any(seed >= mesh.n_points):
        raise ValueError(f"Seed points are outside mesh range 0..{mesh.n_points - 1}")

    config = _read(project_path(defaults.config))
    cfg_lookup = _build_lookup(config)
    n_iter_1, n_iter_2, n_iter_3, top_1, top_2, weight, use_detrain, inverse = load_ev_param(
        crater_name, seed_index, cfg_lookup, cfg_lookup
    )
    # Prefer the exact parameters stored with this APBT result.
    n_iter_1 = int(result.get("N_ITER_1", n_iter_1))
    n_iter_2 = int(result.get("N_ITER_2", n_iter_2))
    n_iter_3 = int(result.get("N_ITER_3", n_iter_3))
    top_1 = float(result.get("TOP_PERCENT_1", top_1))
    top_2 = float(result.get("TOP_PERCENT_2", top_2))
    weight = float(result.get("THINGES_WEIGHT", weight))
    use_detrain = bool(result.get("USE_DETRAIN", use_detrain))
    inverse = bool(result.get("inverse", inverse))
    effective_n_iter_1 = int(result.get("EFFECTIVE_N_ITER_1", n_iter_1))
    effective_weight = float(result.get("EFFECTIVE_THINGES_WEIGHT", weight))
    points_circle = int(result.get("POINTS_CIRCLE", 180))
    ray_length = int(result.get("RAY_LENGTH", 200))

    normals = compute_reference_normals(mesh)
    report("compute_reference_normals -> normals", normals)
    mesh, normals = rotate_mesh_and_normals(mesh, normals, inverse=bool(result.get("inverse", inverse)))
    report("rotate_mesh_and_normals -> mesh", mesh)
    report("rotate_mesh_and_normals -> normals", normals)
    adjacency = adj_vertex_list(mesh)
    report("adj_vertex_list -> adjacency", adjacency)
    working = detrain_mesh(mesh) if use_detrain else mesh
    report("detrain_mesh -> working", working)
    matrix = _based_diffusion_matrix(working, adjacency)
    report("_based_diffusion_matrix -> matrix", matrix)
    best_center = estimate_the_best_center(
        working,
        normals,
        matrix,
        top_1,
        N_ITER_1=effective_n_iter_1,
        THINGES_WEIGHT=effective_weight,
    )
    report("estimate_the_best_center -> best_center", best_center)
    working = return_poly_rim(working, normals, matrix, best_center, smooth_iter=n_iter_2, top_percent_2=top_2)
    report("return_poly_rim -> working", working)
    working, score_field = apply_gradient_diffusion(working, matrix, n_iter=n_iter_3)
    report("apply_gradient_diffusion -> working", working)
    report("apply_gradient_diffusion -> score_field", score_field)
    report("gradient_projection_normalized_smooth", working.point_data[score_field])
    analysis = working.copy()
    best_seeds = find_the_best_seed_points(analysis, seed, score_field, adjacency, n_rings=0)
    report("find_the_best_seed_points -> best_seeds", best_seeds)
    center, radius, shape = find_initial_center(analysis, best_seeds)
    report("find_initial_center -> (center, radius, shape)", (center, radius, shape))
    low_points = find_inside_pts(analysis, center, radius, shape)
    report("find_inside_pts -> low_points", low_points)
    segment_rim = grow_segments_unitl_anchor(analysis, low_points, best_seeds, adjacency)
    report("grow_segments_unitl_anchor -> segment_rim", segment_rim)
    center, radius, rim_shape = _find_center_radius_from_segement_rim(analysis, segment_rim)
    report("_find_center_radius_from_segement_rim -> (center, radius, rim_shape)", (center, radius, rim_shape))
    thinned = get_best_rim_pts_per_bin(analysis, segment_rim, center, score_field, shape_info=rim_shape)
    report("get_best_rim_pts_per_bin -> thinned", thinned)
    smooth = smooth_rim_angular(analysis, thinned, center, score_field, shape_info=rim_shape)
    report("smooth_rim_angular -> smooth", smooth)
    smooth_center = center.copy()
    smooth_center[2] = smooth[:, 2].mean()
    dense = upsample_rim_by_angle(smooth, smooth_center, n_points=points_circle, shape_info=rim_shape)
    report("upsample_rim_by_angle -> dense", dense)
    final_rim = project_polyline_with_rays(working.copy(), segment_rim, dense, radius, ray_length=ray_length)
    report("project_polyline_with_rays -> final_rim", final_rim)
    if len(final_rim) < 2:
        raise ValueError("Could not project enough points for the detected rim")
    return working, final_rim, seed, result, score_field


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize one APBT crater rim.")
    parser.add_argument("--body", required=True, choices=BODIES)
    parser.add_argument("--crater", required=True, type=int, help="Crater number, for example 21")
    parser.add_argument("--seed", type=int, default=0, help="Seed group index, default: 0")
    parser.add_argument("--result", type=Path, help="Optional APBT JSON path")
    parser.add_argument("--save", type=Path, help="Optional screenshot path")
    parser.add_argument("--quiet", action="store_true", help="Hide per-function diagnostic output")
    args = parser.parse_args()

    crater_name = f"crater_{args.crater}"
    mesh, rim, seed, result, score_field = calculate_rim(
        args.body, crater_name, args.seed, args.result, debug=not args.quiet
    )
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, scalars=score_field, cmap="viridis", opacity=0.72, show_scalar_bar=True)
    plotter.add_mesh(_line(rim), color="red", line_width=5, label="APBT detected rim")
    seed_cloud = pv.PolyData(mesh.points[seed])
    plotter.add_mesh(seed_cloud, color="yellow", point_size=14, render_points_as_spheres=True)
    plotter.add_text(
        f"{args.body} {crater_name}  seed={args.seed}\n"
        f"Dmax={result.get('Diameter_max', float('nan')):.2f} m  "
        f"depth={result.get('Depth_med', float('nan')):.2f} m",
        position="upper_left",
        font_size=12,
    )
    plotter.add_legend()
    plotter.show(screenshot=str(project_path(args.save)) if args.save else None)


if __name__ == "__main__":
    main()
