import json
from pathlib import Path


def read_literature(literature_results: str = "inputs/Itokawa/literature.json"):
    base = Path(__file__).resolve().parents[1]
    requested = Path(literature_results)
    candidates = [
        requested,
        base / requested,
        base / "inputs" / "Itokawa" / "literature.json",
    ]

    literature_path = next((path for path in candidates if path.exists()), None)
    if literature_path is None:
        raise FileNotFoundError(
            f"Could not find literature file. Tried: {[str(p) for p in candidates]}"
        )

    with open(literature_path, "r") as f:
        literature_df = json.load(f)

    return literature_df


def _iter_craters(literature_data):
    if isinstance(literature_data, list):
        return literature_data
    if isinstance(literature_data, dict):
        if "crater_candidates" in literature_data:
            return literature_data["crater_candidates"]
        if "crater_database" in literature_data:
            return literature_data["crater_database"]
    raise ValueError(
        "Unsupported literature format. Expected list or dict with "
        "'crater_candidates' / 'crater_database'."
    )


def return_lit_id_naru_hirata():
    literature_df = read_literature("inputs/Itokawa/literature.json")
    literature_by_id_naru_hirata: dict[int, dict] = {}

    for crater in _iter_craters(literature_df):
        crater_id = crater.get("id", crater.get("crater_id"))
        if crater_id is None:
            continue

        diameter = crater.get("diameter_m", {})
        minor_axis = diameter.get("minor_axis_b")
        major_axis = diameter.get("major_axis_a")
        if minor_axis is None or major_axis is None:
            continue

        # Use the central diameter from the reported major/minor axes.
        diameter_mean = float((major_axis + minor_axis) / 2.0)
        diameter_std = float(abs(major_axis - minor_axis) / 2.0)

        literature_by_id_naru_hirata[int(crater_id)] = {
            "naru_hirata_diameter_m": diameter_mean,
            "naru_hirata_diameter_std_m": diameter_std,
            "naru_hirata_minor_axis_b_m": float(minor_axis),
            "naru_hirata_major_axis_a_m": float(major_axis),
            "naru_hirata_depth_m": crater.get("depth_m"),
            "naru_hirata_depth_to_diameter_ratio": crater.get("depth_to_diameter_ratio"),
            "classification": crater.get("confidence_classification"),
            "name": crater.get("name"),
        }

    return literature_by_id_naru_hirata
