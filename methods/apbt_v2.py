import numpy as np
import pyvista as pv
import scipy.sparse as sp

from scipy.spatial.distance import cdist
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression



def _based_diffusion_matrix(poly: pv.PolyData, adj_list: np.ndarray):

    n = poly.n_points
    row, colum, data = [], [], []
    include_self = 1

    for i, neig in enumerate(adj_list):
        neigh_size = len(neig) # normally called deg form degree
        weight = 1.0 / ( neigh_size + include_self )

        row.append(i); colum.append(i); data.append(weight)

        for j in neig:
            row.append(i); colum.append(j); data.append(weight)

    return sp.coo_matrix((data, (row, colum)), shape = (n, n)).tocsr()

def _compute_radial_norm_vector(center: np.ndarray, pts: np.ndarray):
    look_center = 0
    pts_xy = pts[:, :2]
    radial_vector = pts_xy - center[None, :]
    radial_vector_distance = np.linalg.norm(radial_vector, axis = 1)
    radial_vector_distance[radial_vector_distance == look_center] = 1.0
    radial_vector_norm_xy = radial_vector / radial_vector_distance[:, None]
    radial_vector_norm_xyz = np.zeros((pts.shape[0], 3), dtype=float)
    radial_vector_norm_xyz[:, 0:2] = radial_vector_norm_xy

    return radial_vector_norm_xyz

def _projection_smooth(radial_vector_norm_xyz: np.ndarray, normals: np.ndarray, smooth_itr: int = 60, matrix: sp.csr_matrix = None):
    projection = np.sum(radial_vector_norm_xyz * normals, axis = 1)
    projection_smooth = projection.copy()

    for _ in range(smooth_itr):
        projection_smooth = matrix.dot(projection_smooth)

    return projection_smooth

def _compute_gradient_projection(radial_vector_norm_xyz: np.ndarray, grdeint_data: np.ndarray):
    gradient_projection = np.sum(grdeint_data * radial_vector_norm_xyz, axis = 1)
    gradient_projection = np.maximum(gradient_projection, 0.0)
    # gradient_projection = np.abs(gradient_projection)

    return gradient_projection

def _compute_rim_score(poly: pv.PolyData,
                      center_xy: np.ndarray,
                      matrix: sp.csr_matrix,
                      normals: np.ndarray,
                      N_ITER_1: int = 60,
                      ):

    #### VARIABLES for me ####
    temp_name = "rim_projection_smooth"
    smooth_itr = N_ITER_1
    ###
    pts = poly.points
    radial_vector_norm_xyz = _compute_radial_norm_vector(center_xy, pts)
    projection_smooth = _projection_smooth(radial_vector_norm_xyz, normals, smooth_itr, matrix)

    poly.point_data[temp_name] = projection_smooth

    grad_ds = poly.compute_derivative(scalars=temp_name, gradient=True)

    grdeint_data = grad_ds.point_data["gradient"]

    try:
        del poly.point_data[temp_name]

    except Exception:
        print("We can not deltete the name")
        pass

    gradient_projection = _compute_gradient_projection(radial_vector_norm_xyz, grdeint_data)

    return gradient_projection

def _evaluate_ring_quality(poly: pv.PolyData,
                          center_xy: np.ndarray,
                          score: np.ndarray,
                          THINGES_WEIGHT: float = 0.5
                          ):

    top_percent = 97.0
    min_points = 20
    angle_partitions = 36
    thighness_weigth= THINGES_WEIGHT

    pts_xy = poly.points[:, :2]
    threshold = np.percentile(score, top_percent)
    top_pts = np.where(score >= threshold)[0]

    if top_pts.shape[0] < min_points:
        return -np.inf, {"rim_pts": top_pts.size, "coverage": 0.0, "tightness": 0.0}

    vector_to_center = pts_xy[top_pts] - center_xy[None, :]
    distance_to_center = np.linalg.norm(vector_to_center, axis = 1)
    angles = np.arctan2(vector_to_center[:, 1], vector_to_center[:, 0])
    angles_2pi = (angles + 2 * np.pi) % (2 * np.pi)

    angle_index_bins = np.floor(angles_2pi / (2 * np.pi) * angle_partitions).astype(int)
    coverage = np.unique(angle_index_bins).size / angle_partitions

    mean_distance = np.mean(distance_to_center)
    std_distance = np.std(distance_to_center)
    thighness = std_distance / (mean_distance + 1e-9)

    J = coverage - (thighness * thighness_weigth )

    return J,  {"rim_pts": top_pts.size, "coverage": coverage, "tightness": thighness, "thr": threshold}

def _find_first_center(poly: pv.PolyData):
    z = "elev_z"
    low_points: int = 300

    elevation = np.asarray(poly.point_data[z])
    low_points = int(min(max(low_points, 10), poly.n_points))
    low_indices = np.argpartition(elevation, low_points-1)[:low_points]
    center = poly.points[low_indices].mean(axis=0)

    return center

def estimate_the_best_center(poly: pv.PolyData,
                             normals: np.ndarray,
                             diffusion_matrix: sp.csr_matrix,
                             TOP_PERCENT_1: int,
                             N_ITER_1: int = 60,
                             THINGES_WEIGHT: float = 0.5,

                             ):

    max_to_look = TOP_PERCENT_1
    center_factor_search: float = 0.5
    number_of_centers_to_look: int = 10

    center = _find_first_center(poly)
    center_xy = center[:2]
    points_xy = poly.points[:, :2]
    distance_from_center = np.linalg.norm(points_xy - center_xy[None, :], axis = 1)
    max_distance = np.percentile(distance_from_center, max_to_look)
    where_to_look = max_distance * center_factor_search

    if number_of_centers_to_look < 3:
        number_of_centers_to_look = 3

    number_of_centers_to_look = np.linspace(-where_to_look, where_to_look, number_of_centers_to_look)
    candidates = []
    for dx in number_of_centers_to_look:
        for dy in number_of_centers_to_look:
            candidates.append([center_xy[0] + dx, center_xy[1] + dy])

    candidates = np.array(candidates)

    best_center = center_xy.copy()
    best_j = -np.inf

    for center_to_look in candidates:
        score = _compute_rim_score(
            poly,
            center_to_look,
            diffusion_matrix,
            normals,
            N_ITER_1=N_ITER_1,
        )
        J, _ = _evaluate_ring_quality(
            poly,
            center_to_look,
            score,
            THINGES_WEIGHT=THINGES_WEIGHT,
        )

        if J > best_j:
            best_j = J
            best_center = center_to_look

    return best_center

def detrain_mesh(poly: pv.PolyData):
    degree = 2
    points = poly.points
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    xy_stack = np.stack((x, y), axis=1)
    polynomial_base = PolynomialFeatures(degree = degree)
    polynomial_result = polynomial_base.fit_transform(xy_stack)

    model = LinearRegression()
    model.fit(polynomial_result, z)
    print(model.score(polynomial_result, z))
    z_trained= model.predict(polynomial_result)
    z_detrained = z - z_trained

    poly.point_data["z_detrained"] = z_detrained
    poly.point_data["z_trained"] = z_trained

    poly.points[:, 2] = z_detrained
    return poly

def _normalize_projection(projection: np.ndarray, TOP_PERCENT_2: int = 98):
    min_proj = 2

    proj_min, proj_max = np.percentile(projection, [min_proj, TOP_PERCENT_2])
    if proj_max > proj_min:
        projection_normalized = (projection - proj_min) / (proj_max - proj_min)

    projection_normalized = np.clip(projection_normalized, 0.0, 1.0)

    return projection_normalized

def return_poly_rim(poly: pv.PolyData,
                    normals: np.ndarray,
                    matrix: sp.csr_matrix,
                    center_xy: np.ndarray,
                    smooth_iter: int = 60,
                    top_percent_2: int = 98):

    pts = poly.points
    compute_radial_norm_vector_xyz =_compute_radial_norm_vector(center_xy, pts)
    projection_smooth = _projection_smooth(compute_radial_norm_vector_xyz, normals, smooth_iter, matrix)
    temp_name = "rim_projection_smooth"
    poly.point_data[temp_name] = projection_smooth

    gradient = poly.compute_derivative(scalars=temp_name, gradient=True)
    gradient_data = gradient.point_data["gradient"]
    gradient_projection = _compute_gradient_projection(compute_radial_norm_vector_xyz, gradient_data)

    poly.point_data["gradient_projection"] = gradient_projection
    gradient_projection_normalized = _normalize_projection(gradient_projection, top_percent_2)

    poly.point_data["gradient_projection_normalized"] = gradient_projection_normalized

    return poly

def apply_gradient_diffusion(poly: pv.PolyData,
                              matrix: sp.csr_matrix,
                              n_iter: int = 60):

    scalar_base_name = "gradient_projection_normalized"
    scalar_final_name = f"{scalar_base_name}_smooth"

    gradient_projection_normalized = poly.point_data[scalar_base_name]
    gradient_projection_normalized_smooth = gradient_projection_normalized.copy()

    for _ in range(n_iter):
        gradient_projection_normalized_smooth = matrix.dot(gradient_projection_normalized_smooth)

    poly.point_data[scalar_final_name] = gradient_projection_normalized_smooth

    return poly, scalar_final_name

def calculate_rim_diameter(final_rim_pts: np.ndarray):
    min_x = np.min(final_rim_pts[:, 0]); max_x = np.max(final_rim_pts[:, 0])
    min_y = np.min(final_rim_pts[:, 1]); max_y = np.max(final_rim_pts[:, 1])
    diameter_ew = max_x - min_x
    diameter_ns = max_y - min_y
    all_dists = cdist(final_rim_pts, final_rim_pts)
    diameter_max = np.max(all_dists)
    return diameter_ew, diameter_ns, diameter_max






def find_the_best_seed_points(poly: pv.PolyData,
                              seed_points: list,
                              score_field: str,
                              adj_vertexs: np.ndarray,
                              n_rings: int = 0):

    scores = poly.point_data[score_field]

    best = []

    for seed in seed_points:
        neighbors = {seed}
        froniter = {seed}
        for _ in range(n_rings):
            next_frontier = set()

            for pt in froniter:
                for neighbor in adj_vertexs[pt]:
                    if neighbor not in neighbors:

                        neighbors.add(neighbor)
                        next_frontier.add(neighbor)
            froniter = next_frontier

        neighbors = np.array(list(neighbors))
        best_seed = neighbors[np.argmax(scores[neighbors])]
        best.append(best_seed)

    return np.array(best)

def _orthogonal_unit_vector(vec: np.ndarray):
    ortho = np.array([-vec[1], vec[0]], dtype=float)
    norm = np.linalg.norm(ortho)
    if norm < 1e-12:
        return np.array([0.0, 1.0], dtype=float)
    return ortho / norm

def _estimate_shape_from_points_xy(points_xy: np.ndarray, center_xy: np.ndarray = None):
    eps = 1e-9
    points_xy = np.asarray(points_xy, dtype=float)

    if center_xy is None:
        center_xy = points_xy.mean(axis=0)
    else:
        center_xy = np.asarray(center_xy, dtype=float)

    centered = points_xy - center_xy[None, :]

    if centered.shape[0] < 2:
        major_axis = np.array([1.0, 0.0], dtype=float)
        minor_axis = np.array([0.0, 1.0], dtype=float)
        semi_major = 1.0
        semi_minor = 1.0
    else:
        cov = np.cov(centered.T)
        cov = np.atleast_2d(cov)
        if cov.shape != (2, 2):
            cov = np.eye(2, dtype=float)

        eigvals, eigvecs = np.linalg.eigh(cov)
        order = np.argsort(eigvals)[::-1]
        eigvals = np.maximum(eigvals[order], eps)
        eigvecs = eigvecs[:, order]

        major_axis = eigvecs[:, 0]
        major_axis /= (np.linalg.norm(major_axis) + eps)
        minor_axis = _orthogonal_unit_vector(major_axis)

        proj_major = centered @ major_axis
        proj_minor = centered @ minor_axis

        semi_major = max(np.max(np.abs(proj_major)), np.sqrt(eigvals[0]), eps)
        semi_minor = max(np.max(np.abs(proj_minor)), np.sqrt(eigvals[1]), eps)

        if semi_minor > semi_major:
            semi_major, semi_minor = semi_minor, semi_major
            major_axis, minor_axis = minor_axis, major_axis

    axis_ratio = semi_major / max(semi_minor, eps)
    roundness = np.clip(semi_minor / max(semi_major, eps), 0.0, 1.0)
    ellipticity = np.sqrt(max(0.0, 1.0 - (roundness ** 2)))

    return {
        "center_xy": center_xy,
        "major_axis": major_axis,
        "minor_axis": minor_axis,
        "semi_major": semi_major,
        "semi_minor": semi_minor,
        "axis_ratio": axis_ratio,
        "roundness": roundness,
        "ellipticity": ellipticity,
        "is_elliptical": axis_ratio > 1.15,
    }

def _ellipse_normalized_coordinates(points_xy: np.ndarray,
                                    center_xy: np.ndarray,
                                    shape_info: dict):
    eps = 1e-9
    rel = points_xy - center_xy[None, :]
    major_axis = shape_info["major_axis"]
    minor_axis = shape_info["minor_axis"]
    semi_major = max(shape_info["semi_major"], eps)
    semi_minor = max(shape_info["semi_minor"], eps)

    u = rel @ major_axis
    v = rel @ minor_axis

    theta = np.arctan2(v / semi_minor, u / semi_major)
    theta = np.mod(theta, 2.0 * np.pi)
    rho = np.sqrt((u / semi_major) ** 2 + (v / semi_minor) ** 2)
    circular_radius = np.linalg.norm(rel, axis=1)

    return theta, rho, u, v, circular_radius

def _ellipse_xy_from_theta_rho(theta: np.ndarray,
                               rho: np.ndarray,
                               center_xy: np.ndarray,
                               shape_info: dict):
    semi_major = shape_info["semi_major"]
    semi_minor = shape_info["semi_minor"]
    major_axis = shape_info["major_axis"]
    minor_axis = shape_info["minor_axis"]

    local_major = (semi_major * rho * np.cos(theta))[:, None]
    local_minor = (semi_minor * rho * np.sin(theta))[:, None]
    xy = center_xy[None, :] + (local_major * major_axis[None, :]) + (local_minor * minor_axis[None, :])

    return xy

def find_initial_center(poly: pv.PolyData, seed_points: np.ndarray):
    pts_rims = poly.points[seed_points]
    center = pts_rims.mean(axis=0)
    shape_info = _estimate_shape_from_points_xy(pts_rims[:, :2], center[:2])
    center[:2] = shape_info["center_xy"]
    radius = 0.5 * (shape_info["semi_major"] + shape_info["semi_minor"])

    return center, radius, shape_info

def find_inside_pts(poly: pv.PolyData, center: np.ndarray, radius: float, shape_info: dict = None):
    fraq = 0.5
    pts = poly.points
    ptx_xy = pts[:, :2]
    center_xy = center[:2]

    if shape_info is None:
        radius = radius * fraq
        vector_to_center_xy = ptx_xy - center_xy
        distances = np.linalg.norm(vector_to_center_xy, axis=1)
        mask_inside = distances <= radius
    else:
        _, rho, _, _, _ = _ellipse_normalized_coordinates(ptx_xy, center_xy, shape_info)
        mask_inside = rho <= fraq

    inside_pts_indeces = np.where(mask_inside)[0]

    return inside_pts_indeces

def _find_first_rim_points(low_pts: np.ndarray, adj_vertexs: np.ndarray):
    rim_pts = set()
    low_pts_set = set(low_pts)  # Convert to set for faster lookup

    for pt in low_pts:
        for neighbor in adj_vertexs[pt]:
            if neighbor not in low_pts_set:
                rim_pts.add(neighbor)

    return np.array(list(rim_pts))

def _find_closted_low_pt_to_seed_points(poly: pv.PolyData, index_rim_low_pts: np.ndarray, seed_points: np.ndarray):
    low_pts_coordinates = poly.points[index_rim_low_pts]
    seed_points_coordinates = poly.points[seed_points]
    distances = cdist(low_pts_coordinates, seed_points_coordinates)
    closest_seed_points = np.argmin(distances, axis=1)
    return closest_seed_points

def _assign_pts_to_segments(
                            closest_seed_points: np.ndarray,
                            index_rim_low_pts: np.ndarray,
                            seed_points: set,
                            visted_pts: set,
                            where_segment_hits: list,
                            current_segement_rim: list,
                            segments_are_finished: list,
                            final_segement_rim: list):


    for i, index_point in enumerate(index_rim_low_pts):
        segment_id = closest_seed_points[i]
        if index_point in seed_points:
            if not segments_are_finished[segment_id]:
                segments_are_finished[segment_id] = True
                where_segment_hits[segment_id] = index_point
                final_segement_rim[segment_id] = {index_point}
        else:
            current_segement_rim[segment_id].add(index_point)

        visted_pts.add(index_point)
    return visted_pts, where_segment_hits, current_segement_rim, segments_are_finished, final_segement_rim

def _initialize_segments(number_of_segments: int, segment_is_finished: list, current_segement_rim: list, final_segement_rim: list):

    for seg_id in range(number_of_segments):
        if not segment_is_finished[seg_id]:
            if len(current_segement_rim[seg_id]) == 0:
                segment_is_finished[seg_id] = True
            else:
                final_segement_rim[seg_id] = current_segement_rim[seg_id].copy()


    return segment_is_finished, final_segement_rim

def _grow_segements(visited_pts: set,
                    where_segment_hits: list,
                    current_segement_rim: list,
                    segments_are_finished: list,
                    final_segement_rim: list,
                    adj_vertexs: np.ndarray,
                    number_of_segments: int,
                    seed_pts: set
                    ):

    while not all(segments_are_finished):
        any_segment_grew = False
        next_segement_rim = [set() for _ in range(number_of_segments)]

        for seg_id in range(number_of_segments):
            if segments_are_finished[seg_id]: continue

            current_froniter = current_segement_rim[seg_id]
            if len(current_froniter) == 0:
                segments_are_finished[seg_id] = True
                continue

            stopped_this_segment = False

            for pt in current_froniter:
                for neighbor in adj_vertexs[pt]:
                    if neighbor in visited_pts: continue
                    visited_pts.add(neighbor)

                    if neighbor in seed_pts:

                        segments_are_finished[seg_id] = True
                        stopped_this_segment = True
                        where_segment_hits[seg_id] = neighbor
                        break

                    next_segement_rim[seg_id].add(neighbor)
                    any_segment_grew = True

                if stopped_this_segment: break

            if stopped_this_segment:
                current_segement_rim[seg_id] = set()

            else:
                new_froniter = next_segement_rim[seg_id]

                if len(new_froniter) > 0:
                    final_segement_rim[seg_id] = new_froniter.copy()
                    current_segement_rim[seg_id] = new_froniter

                else:
                    segments_are_finished[seg_id] = True
                    current_segement_rim[seg_id] = set()
        if not any_segment_grew and not all(segments_are_finished):
            print("Warning: No segments grew in this iteration, but not all segments are finished. There might be isolated segments."); break


    return [
    np.array(sorted(segment), dtype=np.int64)
    for segment in final_segement_rim
    if len(segment) > 0
]

def grow_segments_unitl_anchor(poly: pv.PolyData, low_pts: np.ndarray, seed_points: np.ndarray, adj_vertexs: np.ndarray):

    final_segments = []
    visited_pts = set(low_pts)
    if len(seed_points) < 4:
        print("We need at least 4 seed points to grow segments until anchor.")
        return []

    number_of_segments = len(seed_points)
    segments_are_finished = [False] * number_of_segments
    where_segment_hits = [None] * number_of_segments
    current_segement_rim = [set() for _ in range(number_of_segments)]
    final_segement_rim = [set() for _ in range(number_of_segments)]
    index_rim_low_pts = _find_first_rim_points(low_pts, adj_vertexs)
    closest_seed_points = _find_closted_low_pt_to_seed_points(poly,
                                                              index_rim_low_pts,
                                                              seed_points)

    (
        visited_pts,
        where_segment_hits,
        current_segement_rim,
        segments_are_finished,
        final_segement_rim,
    ) = _assign_pts_to_segments(
        closest_seed_points,
        index_rim_low_pts,
        set(seed_points),
        visited_pts,
        where_segment_hits,
        current_segement_rim,
        segments_are_finished,
        final_segement_rim,
    )

    segments_are_finished, final_segement_rim = _initialize_segments(number_of_segments,
                                                                     segments_are_finished,
                                                                     current_segement_rim,
                                                                     final_segement_rim)


    final_segments = _grow_segements(
        visited_pts,
        where_segment_hits,
        current_segement_rim,
        segments_are_finished,
        final_segement_rim,
        adj_vertexs,
        number_of_segments,
        set(seed_points)
    )

    if len(final_segments) == 0:
        return np.empty(0, dtype=np.int64)

    return np.unique(np.concatenate(final_segments)).astype(np.int64)

def _find_center_radius_from_segement_rim(poly: pv.PolyData, segment_rim: np.ndarray):
    pts = poly.points[segment_rim]
    center = pts.mean(axis=0)
    shape_info = _estimate_shape_from_points_xy(pts[:, :2], center[:2])
    center[:2] = shape_info["center_xy"]
    radius = 0.5 * (shape_info["semi_major"] + shape_info["semi_minor"])

    return center, radius, shape_info

def get_best_rim_pts_per_bin(poly: pv.PolyData,
                             rim_indices: np.ndarray,
                             center: np.ndarray,
                             score_field: str,
                             n_bins: int = 180,
                             shape_info: dict = None):

    pts = poly.points[rim_indices]
    pts_xy = pts[:, :2]
    center_xy = center[:2]

    if shape_info is None:
        shape_info = _estimate_shape_from_points_xy(pts_xy, center_xy)

    roundness = shape_info["roundness"]
    weight_radial_penalty = 0.06 + (0.06 * roundness)
    n_bins = int(np.clip(round(n_bins * (0.70 + 0.45 * roundness)), 72, 360))

    scores = poly.point_data[score_field][rim_indices]
    theta, rho, _, _, circular_radius = _ellipse_normalized_coordinates(pts_xy, center_xy, shape_info)

    bin_edges = np.linspace(0, 2.0 * np.pi, n_bins + 1)
    best_indices = []

    for bin in range(n_bins):
        is_in_bin = (theta >= bin_edges[bin]) & (theta < bin_edges[bin + 1])
        if not np.any(is_in_bin):
            continue

        local_scores = scores[is_in_bin]
        local_rho = rho[is_in_bin]
        target_rho = np.median(local_rho)
        ellipse_penalty = np.abs(local_rho - target_rho) / (target_rho + 1e-9)

        local_radius = circular_radius[is_in_bin]
        target_radius = np.median(local_radius)
        circular_penalty = np.abs(local_radius - target_radius) / (target_radius + 1e-9)

        combined_penalty = ((1.0 - roundness) * ellipse_penalty) + (roundness * circular_penalty)
        combined_score = local_scores - (weight_radial_penalty * combined_penalty)
        best_local_index = np.argmax(combined_score)

        global_index = np.where(is_in_bin)[0][best_local_index]
        best_indices.append(rim_indices[global_index])

    return np.array(best_indices)

def _sort_values_by_angle(pts: np.ndarray, theta: np.ndarray, radius: np.ndarray, weights: np.ndarray):

    order = np.argsort(theta)
    theta = theta[order]
    radius = radius[order]
    z = pts[:, 2][order]
    weights = weights[order]

    return theta, radius, z, weights

def _angle_weight_calculation(theta: np.ndarray, radius: float, sigma_deg: float = 20.0, weights: np.ndarray = None, numebr_of_iterations: int = 5):

    angle_distance = np.abs(theta[:, None] - theta[None, :])
    angle_best_distance = np.minimum(angle_distance, 2.0 * np.pi - angle_distance)
    distance_weight = np.exp(-(angle_best_distance ** 2) / (2.0 * sigma_deg ** 2))
    final_weight = (distance_weight * (weights[None, :] ** 3))

    for _ in range(numebr_of_iterations):
        radius = (final_weight * radius[None, :]).sum(axis=1) / final_weight.sum(axis=1)

    return radius

def smooth_rim_angular(poly: pv.PolyData,
                       rim_indices: np.ndarray,
                       center: np.ndarray,
                       score_field: str,
                       min_angle: float = 20.0,
                       n_iter: int = 3,
                       shape_info: dict = None):

    clip_min = 0.1
    clip_max = 1.0

    pts = poly.points[rim_indices]
    if shape_info is None:
        shape_info = _estimate_shape_from_points_xy(pts[:, :2], center[:2])

    scores = poly.point_data[score_field][rim_indices]
    s_min, s_max = np.percentile(scores, [3, 98])
    scores_cliped = np.clip((scores - s_min) / (s_max - s_min + 1e-9), clip_min, clip_max)

    angles_2pi, ellipse_radius, _, _, _ = _ellipse_normalized_coordinates(pts[:, :2], center[:2], shape_info)
    angles_2pi, ellipse_radius, z, scores_cliped = _sort_values_by_angle(pts, angles_2pi, ellipse_radius, scores_cliped)

    adaptive_min_angle = np.deg2rad(min_angle * (0.75 + 0.50 * shape_info["roundness"]))
    adaptive_n_iter = max(1, int(round(n_iter * (0.75 + 0.50 * shape_info["roundness"]))))

    ellipse_radius = _angle_weight_calculation(angles_2pi, ellipse_radius.copy(), adaptive_min_angle, scores_cliped, adaptive_n_iter)
    xy_s = _ellipse_xy_from_theta_rho(angles_2pi, ellipse_radius, center[:2], shape_info)
    z_s = z

    return np.column_stack([xy_s[:, 0], xy_s[:, 1], z_s])

def upsample_rim_by_angle(smooth_pts: np.ndarray,
                         smooth_rim_center: np.ndarray,
                         n_points: int = 180,
                         shape_info: dict = None):

    n_points = max(n_points, len(smooth_pts))
    center_xy = smooth_rim_center[:2]

    if shape_info is None:
        shape_info = _estimate_shape_from_points_xy(smooth_pts[:, :2], center_xy)

    theta, ellipse_radius, _, _, _ = _ellipse_normalized_coordinates(smooth_pts[:, :2], center_xy, shape_info)
    z = smooth_pts[:, 2]

    order = np.argsort(theta)
    theta = theta[order]
    ellipse_radius = ellipse_radius[order]
    z = z[order]

    theta_unique, idx_unique = np.unique(theta, return_index=True)
    ellipse_radius = ellipse_radius[idx_unique]
    z = z[idx_unique]
    theta = theta_unique

    theta_ext = np.concatenate([theta, theta[0:1] + 2.0 * np.pi])
    ellipse_radius_ext = np.concatenate([ellipse_radius, ellipse_radius[0:1]])
    z_ext = np.concatenate([z, z[0:1]])

    theta_new = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    ellipse_radius_new = np.interp(theta_new, theta_ext, ellipse_radius_ext)
    z_new = np.interp(theta_new, theta_ext, z_ext)

    xy_new = _ellipse_xy_from_theta_rho(theta_new, ellipse_radius_new, center_xy, shape_info)
    upsmaple_pts = np.column_stack([xy_new[:, 0], xy_new[:, 1], z_new])

    return upsmaple_pts

def project_polyline_with_rays(poly: pv.PolyData, all_rim_indices: np.ndarray,  dense_smooth_pts: np.ndarray, radius: np.ndarray, ray_length: int = 180 ):

    rim_pts_orig = poly.points[all_rim_indices]
    rim_z_ref = np.median(rim_pts_orig[:, 2])

    dense_smooth_for_rays = dense_smooth_pts.copy()
    dense_smooth_for_rays[:, 2] = rim_z_ref

    direction = np.array([0, 0, 1.0])
    ray_length *= radius

    direction /= np.linalg.norm(direction)

    new_pts = []

    for p in dense_smooth_for_rays:
        start = p + direction * ray_length
        end = p - direction * ray_length

        hits, _ = poly.ray_trace(start, end)
        if len(hits):
            new_pts.append(hits[0])

    new_pts = np.asarray(new_pts)

    return new_pts

def frac_neighbors_from_radius(radius: float,
                               r_min: float = 4.0,
                               r_max: float = 200.0,
                               f_min: float = 0.08,
                               f_max: float = 0.50) -> float:

    r = np.clip(radius,r_min,r_max )
    t = (r - r_min) / (r_max - r_min)
    frac = f_min + t * (f_max - f_min)
    #np.interp(radius [r_min, r_max],[f_min, f_max])
    return frac
