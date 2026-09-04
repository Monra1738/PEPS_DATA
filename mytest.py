"""Save top-view APBT projection fields for crater meshes."""

import argparse
import json

import numpy as np
import pyvista as pv

from methods.apbt_v2 import (
    _based_diffusion_matrix,
    _find_center_radius_from_segement_rim,
    apply_gradient_diffusion,
    detrain_mesh,
    estimate_the_best_center,
    find_initial_center,
    find_inside_pts,
    find_the_best_seed_points,
    get_best_rim_pts_per_bin,
    grow_segments_unitl_anchor,
    return_poly_rim,
    smooth_rim_angular,
)
from scripts.batch_config import APBT_PATHS, project_path
from utils.mesh_manipulation import (
    adj_vertex_list,
    compute_reference_normals,
    force_meters,
    rotate_mesh_and_normals,
)


def calculate_projection(body, crater, paths, config, seed_indices):
    selected_cfg = config[f"{body}/{crater}"][-1]
    n_iter_1 = int(selected_cfg["n_iter_1"])
    n_iter_2 = int(selected_cfg["n_iter_2"])
    n_iter_3 = int(selected_cfg["n_iter_3"])
    top_percent_1 = float(selected_cfg["top_percent_1"])
    top_percent_2 = float(selected_cfg["top_percent_2"])
    things_weight = float(selected_cfg["things_weight"])
    use_detrain = bool(selected_cfg["USE_DETRAIN"])
    inverse = bool(selected_cfg["INVERSE"])

    mesh_path = project_path(paths.mesh_dir) / f"{crater}.vtk"
    mesh = force_meters(pv.read(mesh_path))
    normals = compute_reference_normals(mesh)
    mesh, normals = rotate_mesh_and_normals(mesh, normals, inverse=inverse)
    adjacency = adj_vertex_list(mesh)

    projection_mesh = detrain_mesh(mesh) if use_detrain else mesh
    matrix = _based_diffusion_matrix(projection_mesh, adjacency)
    best_center = estimate_the_best_center(
        projection_mesh,
        normals,
        matrix,
        top_percent_1,
        N_ITER_1=n_iter_1,
        THINGES_WEIGHT=things_weight,
    )
    projection_mesh = return_poly_rim(
        projection_mesh,
        normals,
        matrix,
        best_center,
        smooth_iter=n_iter_2,
        top_percent_2=top_percent_2,
    )
    projection_mesh, score_field = apply_gradient_diffusion(
        projection_mesh,
        matrix,
        n_iter=n_iter_3,
    )

    best_seeds = find_the_best_seed_points(
        projection_mesh,
        seed_indices,
        score_field,
        adjacency,
        n_rings=0,
    )
    center, radius, seed_shape = find_initial_center(projection_mesh, best_seeds)
    low_points = find_inside_pts(projection_mesh, center, radius, seed_shape)
    segment_rim = grow_segments_unitl_anchor(
        projection_mesh,
        low_points,
        best_seeds,
        adjacency,
    )
    center, _, rim_shape = _find_center_radius_from_segement_rim(
        projection_mesh,
        segment_rim,
    )
    thinned_rim = get_best_rim_pts_per_bin(
        projection_mesh,
        segment_rim,
        center,
        score_field,
        shape_info=rim_shape,
    )
    smooth_points = smooth_rim_angular(
        projection_mesh,
        thinned_rim,
        center,
        score_field,
        shape_info=rim_shape,
    )
    return projection_mesh, score_field, smooth_points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--body", default="all", choices=["all", *APBT_PATHS])
    parser.add_argument("--crater", type=int, help="Only process this crater")
    parser.add_argument("--seed", type=int, default=0, help="Seed group to display")
    parser.add_argument("--output", default="projection_visualizations")
    args = parser.parse_args()

    bodies = APBT_PATHS if args.body == "all" else [args.body]
    output = project_path(args.output)

    for body in bodies:
        paths = APBT_PATHS[body]
        with project_path(paths.config).open() as handle:
            config = json.load(handle)
        with project_path(paths.seeds).open() as handle:
            seeds = json.load(handle)

        meshes = sorted(
            project_path(paths.mesh_dir).glob("crater_*.vtk"),
            key=lambda path: int(path.stem.rsplit("_", 1)[1]),
        )
        for mesh_path in meshes:
            crater = mesh_path.stem
            if args.crater is not None and crater != f"crater_{args.crater}":
                continue
            if f"{body}/{crater}" not in config:
                continue
            seed_groups = seeds.get(crater, [])
            if args.seed < 0 or args.seed >= len(seed_groups):
                print(f"skipped {body} {crater}: seed {args.seed} not found", flush=True)
                continue
            seed_indices = np.asarray(seed_groups[args.seed], dtype=int)

            print(f"processing {body} {crater}", flush=True)
            projection_mesh, score_field, smooth_points = calculate_projection(
                body, crater, paths, config, seed_indices
            )

            plotter = pv.Plotter(off_screen=True, window_size=(1200, 1200))
            plotter.add_mesh(
                projection_mesh,
                scalars=score_field,
                cmap="terrain",
                scalar_bar_args={"title": score_field},
            )
            # plotter.add_mesh(
            #     pv.PolyData(projection_mesh.points[seed_indices]),
            #     color="yellow",
            #     point_size=14,
            #     render_points_as_spheres=True,
            # )
            # smooth_loop = np.vstack((smooth_points, smooth_points[0]))
            # plotter.add_mesh(
            #     pv.lines_from_points(smooth_loop),
            #     color="red",
            #     line_width=5,
            # )
            # plotter.add_mesh(
            #     pv.PolyData(smooth_points),
            #     color="red",
            #     point_size=7,
            #     render_points_as_spheres=True,
            # )
            # plotter.add_text(f"{body} {crater} seed={args.seed}", font_size=12)
            plotter.view_xy()
            plotter.camera.parallel_projection = True
            destination = output / body / f"{crater}_seed_{args.seed}.png"
            destination.parent.mkdir(parents=True, exist_ok=True)
            plotter.show(screenshot=str(destination), auto_close=True)
            plotter.close()
            print(f"saved {destination}", flush=True)


if __name__ == "__main__":
    main()
