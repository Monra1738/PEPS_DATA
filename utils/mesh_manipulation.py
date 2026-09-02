import pyvista as pv
import numpy as np
from sklearn.decomposition import PCA
from matplotlib.path import Path

def rotate_mesh(poly : pv.PolyData, inverse = False):
    vertices = poly.points
    pca = PCA(n_components=3)
    pca.fit(vertices)
    principal_axes = pca.components_
    if principal_axes[2,2]<0:
        rotation_matrix = np.array([principal_axes[0], principal_axes[1], -principal_axes[2]])
    else:
        rotation_matrix = np.array([principal_axes[0], principal_axes[1], principal_axes[2]])
    rotated_vertices = np.dot(vertices, rotation_matrix.T)
    poly.points = rotated_vertices
    pts = poly.points
    if inverse :
        pts[:,2] *= -1
    poly.points = pts
    z_coords = poly.points[:, 2]
    poly.point_data['elev_z'] = z_coords
    return poly

def adj_vertex_list(poly : pv.PolyData):
    n_pts = poly.n_points
    arr_faces = poly.faces
    neigh = [set() for _ in range(n_pts)]
    i,len_face = 0,len(arr_faces)
    while (i< len_face ):
        m = arr_faces[i]; i +=1
        triangle = arr_faces[i:i+m]; i+=m
        try:
            for indx in range(m):
                v = triangle[indx]; u = triangle[(indx+1)%m]
                if v != u:
                    neigh[v].add(u);neigh[u].add(v)
        except:
            return None
    return [list(s) for s in neigh]

def pooling_values(values : np.ndarray, poly : pv.PolyData, rings : int = 2, method = "median" ) -> np.ndarray:

    neigh = adj_vertex_list(poly)
    result = np.empty_like(values)
    n_points = poly.n_points
    for pt_i in range(n_points):
        visited = {pt_i};frontire = {pt_i}
        for _ in range(rings):
            nxt = set()
            for f in frontire: nxt.update(neigh[f])
            nxt -= visited
            visited.update(nxt)
            frontire = nxt
            if not frontire: break
        v = values[list(visited)]
        if method == "median": result[pt_i] = float(np.median(v))
        else: result[pt_i] = float(np.mean(v))
    return result

def signal_choice( H_s : np.ndarray, k1 : np.ndarray,k2 : np.ndarray, s_choice = "H"):
    if s_choice == "H": signal = H_s
    elif  s_choice == "k1": signal = k1
    elif s_choice == "absH": signal = np.abs(H_s)
    elif s_choice == "abs_k1": signal = np.abs(k1)
    elif s_choice == "abs_k2": signal = np.abs(k2)
    return signal

def ensure_normals(poly : pv.PolyData):
    poly.compute_normals(cell_normals=True, point_normals=True, split_vertices= False,
                         consistent_normals=True, auto_orient_normals=True,
                         inplace=True, feature_angle=180)
    pn = poly.point_data.get('Normals', None); cn = poly.cell_data.get('Normals')
    if pn is None: pn = poly.point_normals
    if cn is None: cn = poly.cell_normals
    poly.cell_data['CellNormals'] = np.ascontiguousarray(cn)
    poly.point_data['PointNormals'] = np.ascontiguousarray(pn)
    cn = poly.cell_data.get("Normals", None)
    poly.cell_data['CellNormals'] = np.ascontiguousarray(cn)
    return poly

def curv_copy_mesh(poly : pv.PolyData,taub_n_iter, taubin_p_band,  pool_ring ,pool_method ):
    mesh_copy = poly.copy()
    mesh_copy = mesh_copy.smooth_taubin(n_iter = taub_n_iter , pass_band = taubin_p_band)
    mesh_copy = ensure_normals(mesh_copy)
    H_raw = mesh_copy.curvature("Mean"); K_raw = mesh_copy.curvature("Gaussian")
    H_s = pooling_values(H_raw, mesh_copy, rings=pool_ring, method=pool_method)
    K_s = pooling_values(K_raw, mesh_copy, rings=pool_ring, method=pool_method)
    Delta = np.maximum(H_s**2 - K_s, 0.0); root_D = np.sqrt(Delta)
    k1 = H_s +root_D; k2 = H_s - root_D
    k1, k2 = np.maximum(k1, k2), np.minimum(k1, k2)
    return H_raw, K_raw, H_s, K_s, k1, k2

def force_meters(poly: pv.PolyData, scale_hint = None) -> pv.PolyData:

    s = None; pts = poly.points
    if isinstance(scale_hint, (int, float)): s = float(scale_hint)
    elif scale_hint in (None, "m"): s = 1000
    elif scale_hint == "mm": s = 1e-6
    elif scale_hint == "km": s = 1e6
    else:
        span = np.ptp(pts, axis =0).max()
        s = 1e-6 if span > 1e5 else 1.0
    if s != 1.0:
        pts*=s
        poly.points = pts
    poly.field_data["unit_scale_to_m"] = np.array([s])
    return poly


def scale_from_source_unit_to_meters(
    poly: pv.PolyData,
    source_unit: str,
) -> pv.PolyData:

    normalized_unit = str(source_unit).strip().lower()
    factors = {
        "m": 1.0,
        "meter": 1.0,
        "meters": 1.0,
        "metre": 1.0,
        "metres": 1.0,
        "km": 1000.0,
        "kilometer": 1000.0,
        "kilometers": 1000.0,
        "kilometre": 1000.0,
        "kilometres": 1000.0,
        "cm": 0.01,
        "mm": 0.001,
    }
    if normalized_unit not in factors:
        raise ValueError(
            f"Unsupported source unit {source_unit!r}; "
            f"expected one of {sorted(factors)}"
        )

    factor = factors[normalized_unit]
    if factor != 1.0:
        poly.points = np.asarray(poly.points, dtype=float) * factor
    poly.field_data["unit_scale_to_m"] = np.array([factor])
    poly.field_data["source_unit"] = np.array([normalized_unit])
    return poly

def frac_neighbors_from_radius(radius: float,
                               r_min: float = 4.0,
                               r_max: float = 200.0,
                               f_min: float = 0.08,
                               f_max: float = 0.30) -> float:

    r = np.clip(radius,r_min,r_max )
    t = (r - r_min) / (r_max - r_min)
    frac = f_min + t * (f_max - f_min)
    #np.interp(radius [r_min, r_max],[f_min, f_max])
    return frac

def return_depth(poly: pv.PolyData,
                 rim_pts: np.ndarray,
                 center_xy: np.ndarray,
                 frac_neighbors: float):

    rim_pts   = np.asarray(rim_pts)
    center_xy = np.asarray(center_xy).reshape(2)
    if rim_pts.ndim != 2 or rim_pts.shape[1] != 3:

        raise ValueError(f"rim_pts must be (N,3), got {rim_pts.shape}")

    rim_z = rim_pts[:, 2]
    rim_level = np.median(rim_z)

    rim_xy = rim_pts[:,:2]
    rim_poly_2d = Path(rim_xy)

    pts = poly.points
    pts_xy = pts[:, :2]

    inside_mask = rim_poly_2d.contains_points(pts_xy)
    inside_pts = pts[inside_mask]



    if inside_pts.shape[0] == 0:
        print("[WARN] No points found inside rim polygon.")
        return None, None, None, None, None


    d = np.linalg.norm(inside_pts[:, :2] - center_xy[None, :], axis=1)
    n = max(1, int(frac_neighbors * inside_pts.shape[0]))

    nearest_idx = np.argpartition(d, n-1)[:n]
    nearest_pts = inside_pts[nearest_idx]

    floor_z = nearest_pts[:, 2]
    floor_level = float(np.min(floor_z))

    depth = rim_level - floor_level
    depths_per_rim = rim_z - floor_level

    depth_med = float(np.median(depths_per_rim))
    depth_p40, depth_p60 = np.percentile(depths_per_rim, [40, 60])

    err_minus = depth_med - depth_p40
    err_plus = depth_p60 - depth_med

    # print(f"---- Crater depth (polygon interior, {int(frac_neighbors*100)}% nearest) ----")
    # print(f"Rim level (median z)     : {rim_level:.3f}")
    # print(f"Floor level (min z)      : {floor_level:.3f}")
    # print(f"Depth (median)           : {depth_med:.3f}")
    # print(f"Depth                    : {depth:.3f}")
    # print(f"Depth P40,P60            : {depth_p40:.3f}, {depth_p60:.3f}")
    # print(f"Depth error (+/-)        : -{err_minus:.3f} / +{err_plus:.3f}")
    # print(f"Interior pts total       : {inside_pts.shape[0]}")
    # print(f"Used nearest {n} pts     : (frac={frac_neighbors})")

    return depth_med, depth,depth_p40,depth_p60, rim_level, floor_level, err_minus, err_plus, frac_neighbors

def compute_reference_normals(poly: pv.PolyData):
    poly = poly.copy(deep=True)
    poly.compute_normals(
        cell_normals=False,
        point_normals=True,
        split_vertices=False,
        consistent_normals=True,
        auto_orient_normals=False,   # important
        inplace=True,
        feature_angle=180,
    )
    normals = np.asarray(poly.point_data["Normals"], dtype=float)
    normals /= (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-12)
    return normals

def rotate_mesh_and_normals(poly: pv.PolyData,
                            normals: np.ndarray,
                            inverse: bool = False):
    poly = poly.copy(deep=True)
    vertices = np.asarray(poly.points, dtype=float)

    pca = PCA(n_components=3)
    pca.fit(vertices)

    R = pca.components_.copy()
    if R[2, 2] < 0:
        R[2] *= -1.0

    A = R.copy()

    # this reproduces your "inverse" behavior on points
    if inverse:
        A[2, :] *= -1.0

    rotated_vertices = vertices @ A.T
    rotated_normals = np.asarray(normals, dtype=float) @ A.T

    rotated_normals /= (np.linalg.norm(rotated_normals, axis=1, keepdims=True) + 1e-12)

    poly.points = rotated_vertices
    poly.point_data["Normals"] = rotated_normals.astype(np.float32)
    poly.point_data["elev_z"] = poly.points[:, 2]

    return poly, rotated_normals
