import json
import numpy as np
from pathlib import Path


def read_literature(literature_results: str = "inputs/Bennu/literature.json"):
    base = Path(__file__).resolve().parents[1]
    requested = Path(literature_results)
    candidates = [
        requested,
        base / requested,
        base / "inputs" / "Bennu" / "literature.json",
    ]

    literature_path = next((path for path in candidates if path.exists()), None)
    if literature_path is None:
        raise FileNotFoundError(
            f"Could not find literature file. Tried: {[str(p) for p in candidates]}"
        )

    with open(literature_path, "r") as f:
        literature_v2_df = json.load(f)

    return literature_v2_df


def _iter_craters(literature_data):
    if isinstance(literature_data, list):
        return literature_data
    if isinstance(literature_data, dict) and "crater_database" in literature_data:
        return literature_data["crater_database"]
    raise ValueError("Unsupported literature format. Expected list or dict with crater_database.")


def return_lit_id_deshapriya():
    literature_v2_df = read_literature("inputs/Bennu/literature.json")
    literature_by_id_deshapriya: dict[int, dict] = {}

    for crater in _iter_craters(literature_v2_df):
        deshapriya = crater.get("deshapriya")
        if deshapriya is None:
            continue

        c_id = crater["crater_id"]
        deshapriya_diam = deshapriya.get("diameter_m")
        if deshapriya_diam is None:
            continue

        literature_by_id_deshapriya[c_id] = {
            "deshapriya_diameter_m": deshapriya_diam,
            "deshapriya_diameter_std_m": deshapriya.get("diameter_std_m", 0.0),
        }

    return literature_by_id_deshapriya


def return_lit_id_daly():
    literature_v2_df = read_literature("inputs/Bennu/literature.json")
    literature_by_id_daly: dict[int, dict] = {}

    for crater in _iter_craters(literature_v2_df):
        daly = crater.get("daly")
        if daly is None:
            continue

        c_id = crater["crater_id"]
        daly_diam = daly.get("diameter_m")
        if daly_diam is None:
            continue

        literature_by_id_daly[c_id] = {
            "daly_diameter_mean": daly_diam,
            "daly_diameter_std": daly.get("diameter_std_m"),
            "daly_depth_mean": daly.get("depth_m"),
            "daly_depth_std": daly.get("depth_std_m"),
            "daly_ratio_mean": daly.get("d_over_d"),
            "daly_ratio_std": daly.get("d_over_d_std"),
        }

    return literature_by_id_daly


def return_lit_id_bierhaus():
    literature_v2_df = read_literature("inputs/Bennu/literature.json")
    literature_by_id_bierhaus: dict[int, dict] = {}

    for crater in _iter_craters(literature_v2_df):
        bierhaus = crater.get("bierhaus")
        if bierhaus is None:
            continue

        c_id = crater["crater_id"]
        bierhaus_diam = bierhaus.get("diameter_m")
        if bierhaus_diam is None:
            continue

        literature_by_id_bierhaus[c_id] = {
            "bierhaus_diameter_mean": bierhaus_diam,
        }
    print()

    return literature_by_id_bierhaus


def compare_daly_bierhaus_deshapriya() -> dict:
    daly = return_lit_id_daly()
    bierhaus = return_lit_id_bierhaus()
    deshapriya = return_lit_id_deshapriya()
    result = {}

    for crater_id, crater_info in deshapriya.items():
        if crater_id not in daly or crater_id not in bierhaus:
            continue

        deshapriya_std = crater_info["deshapriya_diameter_std_m"]
        deshapriya_diam_mean = crater_info["deshapriya_diameter_m"]
        deshapriya_min_diam = deshapriya_diam_mean - deshapriya_std
        deshapriya_max_diam = deshapriya_diam_mean + deshapriya_std

        daly_diam = daly[crater_id]["daly_diameter_mean"]
        bierhaus_diam = bierhaus[crater_id]["bierhaus_diameter_mean"]

        mean_err = np.abs(deshapriya_diam_mean - daly_diam)
        gap_error_daly_deshapriya = 0.0
        gap_error_bierhaus_deshapriya = 0.0

        if daly_diam < deshapriya_min_diam:
            gap_error_daly_deshapriya = deshapriya_min_diam - daly_diam
        elif daly_diam > deshapriya_max_diam:
            gap_error_daly_deshapriya = daly_diam - deshapriya_max_diam

        if bierhaus_diam < deshapriya_min_diam:
            gap_error_bierhaus_deshapriya = deshapriya_min_diam - bierhaus_diam
        elif bierhaus_diam > deshapriya_max_diam:
            gap_error_bierhaus_deshapriya = bierhaus_diam - deshapriya_max_diam

        gap_error_bierhaus_daly = np.abs(bierhaus_diam - daly_diam)

        result[crater_id] = {
            "mean_err": mean_err,
            "gap_error_daly_deshapriya": gap_error_daly_deshapriya,
            "gap_error_bierhaus_deshapriya": gap_error_bierhaus_deshapriya,
            "daly_diam": daly_diam,
            "deshapriya_diam": deshapriya_diam_mean,
            "deshapriya_std": deshapriya_std,
            "bierhaus_diam": bierhaus_diam,
            "gap_error_bierhaus_daly": gap_error_bierhaus_daly,
        }

    return result
