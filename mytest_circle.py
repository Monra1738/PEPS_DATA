"""Save top-view CIRCLE fields and projected rims for crater meshes."""

import argparse
import json
from pathlib import Path

import numpy as np
import pyvista as pv

from methods.circle import find_crater_rim_from_seeds, project_circle_with_rays
from scripts.batch_config import CIRCLE_PATHS, project_path
from utils.mesh_manipulation import (
    curv_copy_mesh,
    force_meters,
    rotate_mesh,
    signal_choice,
)


def calculate_circle(body, crater, paths, config, seed_indices, seed_index):
    config_value = next(
        value for key, value in config.items() if Path(key).stem == crater
    )
    configurations = config_value if isinstance(config_value, list) else [config_value]
    selected_cfg = configurations[min(seed_index, len(configurations) - 1)]
    if "TAUBIN_N_ITER" in selected_cfg:
        taubin_iterations = int(selected_cfg.get("TAUBIN_N_ITER", 1))
        taubin_pass_band = float(selected_cfg.get("TAUBIN_PASS_BAND", 0.1))
        pool_rings = int(selected_cfg.get("POOL_RINGS", 1))
        pool_method = str(selected_cfg.get("POOL_METHOD", "median"))
        signal_name = str(selected_cfg.get("SIGNAL_CHOICE", "abs_k1"))
        inverse = bool(selected_cfg.get("INVERSE", False))
    else:
        taubin_iterations = int(selected_cfg["Smoothing: "]["n_iter"])
        taubin_pass_band = float(selected_cfg["Smoothing: "]["pass_band"])
        pool_rings = int(selected_cfg["Pooling : "]["number_rings"])
        pool_method = str(selected_cfg["Pooling : "]["pool_method"])
        signal_name = str(selected_cfg["Signal"])
        inverse = bool(selected_cfg["Mesh_inverse"])

    mesh_path = project_path(paths.mesh_dir) / f"{crater}.vtk"
    mesh = rotate_mesh(pv.read(mesh_path), inverse)
    mesh = force_meters(mesh)
    _, _, mean_curvature, _, k1, k2 = curv_copy_mesh(
        mesh,
        taubin_iterations,
        taubin_pass_band,
        pool_rings,
        pool_method,
    )
    circle_signal = signal_choice(
        mean_curvature,
        k1,
        k2,
        s_choice=signal_name,
    )
    mesh.point_data["circle_signal"] = circle_signal

    center_xy, radius, _, _, rim_indices = find_crater_rim_from_seeds(
        mesh,
        seed_indices,
        circle_signal,
        0.25,
        0.70,
        0.99,
    )
    if center_xy is None:
        raise ValueError("CIRCLE could not find the crater rim")

    rim_points = np.asarray(mesh.points[rim_indices])
    theta = np.linspace(0, 2 * np.pi, 200, endpoint=True)
    circle_points = np.column_stack(
        (
            center_xy[0] + radius * np.cos(theta),
            center_xy[1] + radius * np.sin(theta),
            np.full_like(theta, np.median(rim_points[:, 2])),
        )
    )
    projected_rim = np.asarray(
        project_circle_with_rays(
            mesh,
            circle_points,
            1e3,
            direction=(0, 0, -1),
        )
    )
    return mesh, projected_rim, rim_points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--body", default="all", choices=["all", *CIRCLE_PATHS])
    parser.add_argument("--crater", type=int, help="Only process this crater")
    parser.add_argument("--seed", type=int, default=0, help="Seed group to display")
    parser.add_argument("--output", default="circle_visualizations")
    args = parser.parse_args()

    bodies = CIRCLE_PATHS if args.body == "all" else [args.body]
    output = project_path(args.output)

    for body in bodies:
        paths = CIRCLE_PATHS[body]
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
            if not any(Path(key).stem == crater for key in config):
                continue
            seed_groups = seeds.get(crater, [])
            if args.seed < 0 or args.seed >= len(seed_groups):
                print(f"skipped {body} {crater}: seed {args.seed} not found", flush=True)
                continue
            seed_indices = np.asarray(seed_groups[args.seed], dtype=int)

            print(f"processing {body} {crater} seed={args.seed}", flush=True)
            try:
                mesh, projected_rim, rim_points = calculate_circle(
                    body,
                    crater,
                    paths,
                    config,
                    seed_indices,
                    args.seed,
                )
            except Exception as error:
                print(f"skipped {body} {crater}: {error}", flush=True)
                continue

            plotter = pv.Plotter(off_screen=True, window_size=(1200, 1200))
            plotter.add_mesh(
                mesh,
                scalars="circle_signal",
                cmap="terrain",
                scalar_bar_args={"title": "circle_signal"},
            )
            # plotter.add_mesh(
            #     pv.PolyData(mesh.points[seed_indices]),
            #     color="yellow",
            #     point_size=14,
            #     render_points_as_spheres=True,
            # )
          
            # rim_loop = np.vstack((projected_rim, projected_rim[0]))
            # plotter.add_mesh(
            #     pv.lines_from_points(rim_loop),
            #     color="red",
            #     line_width=5,
            # )
            plotter.view_xy()
            plotter.camera.parallel_projection = True
            plotter.add_text(f"CIRCLE {body} {crater} seed={args.seed}", font_size=12)
            destination = output / body / f"{crater}_seed_{args.seed}.png"
            destination.parent.mkdir(parents=True, exist_ok=True)
            plotter.show(screenshot=str(destination), auto_close=True)
            plotter.close()
            print(f"saved {destination}", flush=True)


if __name__ == "__main__":
    main()
