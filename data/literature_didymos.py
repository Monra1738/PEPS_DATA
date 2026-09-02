import json
from pathlib import Path


def read_literature(literature_results: str = "inputs/Didymos/literature.json"):
    base = Path(__file__).resolve().parents[1]
    requested = Path(literature_results)
    candidates = [
        requested,
        base / requested,
        base / "inputs" / "Didymos" / "literature.json",
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
        if "craters" in literature_data:
            return literature_data["craters"]
        if "crater_candidates" in literature_data:
            return literature_data["crater_candidates"]
        if "crater_database" in literature_data:
            return literature_data["crater_database"]
    raise ValueError(
        "Unsupported literature format. Expected list or dict with "
        "'craters' / 'crater_candidates' / 'crater_database'."
    )


def return_lit_id_barnouin():
    literature_df = read_literature("inputs/Didymos/literature.json")
    literature_by_id_Barnouin: dict[int, dict] = {}

    for crater in _iter_craters(literature_df):
        crater_id = crater.get("crater_id", crater.get("id"))
        if crater_id is None:
            continue

        # Current Didymos file provides one diameter value per crater.
        diameter_m = crater.get("diameter_m")
        if diameter_m is None:
            continue

        location = crater.get("location", {})
        literature_by_id_Barnouin[int(crater_id)] = {
            "Barnouin_diameter_m": float(diameter_m),
            "Barnouin_diameter_std_m": 0.0,
            "classification": crater.get("confidence"),
            "name": crater.get("id"),
            "latitude": location.get("latitude"),
            "longitude": location.get("longitude"),
        }

    return literature_by_id_Barnouin
