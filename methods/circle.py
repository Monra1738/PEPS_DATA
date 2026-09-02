
import pyvista as pv
import numpy as np
import math

def fit_circle(seed_xy: np.ndarray):
    x = seed_xy[:, 0]
    y = seed_xy[:, 1]
    A = np.c_[x, y, np.ones_like(x)]
    b = -(x**2 + y**2)

    D, E, F = np.linalg.lstsq(A, b, rcond=None)[0]
    cx, cy = -D/2.0, -E/2.0
    rad2 = cx**2 + cy**2 - F
    rad = math.sqrt(rad2) if rad2 > 0 else 0.0
    return np.array([cx, cy]), float(rad)

def find_crater_rim_from_seeds(poly: pv.PolyData, seed_ids : list, curvature_values: np.ndarray, band_rim: 0.25,
                               curv_low: 0.7, curv_high: 0.98):
    poly_pts = poly.points
    seed_ids = np.asarray(seed_ids, dtype= int)

    if poly_pts.shape[0] != curvature_values.shape[0]:
        raise ValueError("curvature_values length must match poly.points length")


    seed_xy = poly_pts[seed_ids, :2]
    center_ini_xy, radi_ini = fit_circle(seed_xy)

    poly_points_xy  = poly_pts[:, :2]
    radi_all_pts = np.linalg.norm(poly_points_xy - center_ini_xy[None,:], axis = 1)

    band = radi_ini * band_rim

    candidate_mask = np.abs(radi_all_pts - radi_ini) <= band


    if candidate_mask.sum() < 10: candidate_mask = np.abs(radi_all_pts - radi_ini) <= (band * 2)
    if candidate_mask.sum() < 3: print("Can not find enouth rim candidates. Use the other methdod "); return  None, None, None, None, None

    signal = np.asarray(curvature_values)
    signal_candidate = signal[candidate_mask]
    low = np.quantile(signal_candidate, curv_low)
    high = np.quantile(signal_candidate, curv_high)

    rim_mask = candidate_mask & (curvature_values >= low) & (curvature_values <= high)

    if rim_mask.sum() <= 3: rim_mask = candidate_mask & (curvature_values >= low)

    rim_indices = np.where(rim_mask)[0]

    if rim_indices.size <=3: print("fit_crater_rim_from_seeds: You should change the low andhigh '%' curv values. ");return None, None, None, None, None #size here


    rim_pts_xy = poly_points_xy[rim_indices]

    center_xy, final_radius = fit_circle(rim_pts_xy)
    dist_center = np.linalg.norm(poly_points_xy - center_xy[None, :], axis=1)

    center_id = int(np.argmin(dist_center))
    diamter = final_radius*2

    print(f"Cercle refinat (curvatura): centre≈({center_xy[0]:.2f}, {center_xy[1]:.2f}), R≈{final_radius:.2f}")
    print(f"Node més proper com a nou centre: {center_id}, diàmetre≈{diamter:.2f}")

    return center_xy, final_radius, center_id, diamter, rim_indices

def project_circle_with_rays(poly : pv.PolyData, pts_xyz: np.ndarray, ray_length, direction):

    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)

    new_pts = []
    for p in pts_xyz:
        start = p + direction*ray_length
        end  = p - direction*ray_length
        hits,_ = poly.ray_trace(start,end)
        if len(hits):
            new_pts.append(hits[0])


    new_pts = np.asarray(new_pts) #better for iteration nparray than just a  list
    return new_pts
