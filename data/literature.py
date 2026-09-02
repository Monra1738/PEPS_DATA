import json
import numpy as np
from pathlib import Path


def read_literature(literature_results: str = 'inputs/Ryugu/literature.json'):
    base = Path(__file__).resolve().parents[1]
    requested = Path(literature_results)
    candidates = [
        requested,
        base / requested,
        base / "inputs" / "Ryugu" / "literature.json",
    ]

    literature_path = next((path for path in candidates if path.exists()), None)
    if literature_path is None:
        raise FileNotFoundError(
            f"Could not find literature file. Tried: {[str(p) for p in candidates]}"
        )

    with open(literature_path, "r") as f:
        literature_v2_df = json.load(f)

    return literature_v2_df

def return_lit_id_hirata():

    literature_v2_df = read_literature('inputs/Ryugu/literature.json')
    literature_by_id_hirata: dict[int, dict] = {}

    for crater in literature_v2_df["crater_database"]:
        c_id = crater["crater_id"]

        hirata = crater["hirata_2020"]
        hirata_diam = hirata["diameter_m"]
        hirata_classification = hirata["classification"]

        literature_by_id_hirata[c_id] = {
                "hirata_diameter_m": hirata_diam,
                "classification" : hirata_classification,
            }
    return literature_by_id_hirata

def return_lit_id_noguchi():
    literature_by_id_noguchi: dict[int, dict] = {}
    literature_v2_df = read_literature('inputs/Ryugu/literature.json')
    literature_by_id_hirata: dict[int, dict] = {}

    for crater in literature_v2_df["crater_database"]:
        c_id = crater["crater_id"]

        noguchi = crater.get("noguchi_2021")
        noguchi_diam = noguchi["diameter_m"]
        noguchi_depth = noguchi["depth_m"]
        noguchi_ratio = noguchi["depth_to_diameter_ratio"]
        noguchi_classification = noguchi["classification"]

        literature_by_id_noguchi[c_id] = {
        "noguchi_diameter_mean": noguchi_diam["mean"],
            "noguchi_diameter_std": noguchi_diam["std_dev"],
            "noguchi_depth_mean": noguchi_depth["mean"],
            "noguchi_depth_std": noguchi_depth["std_dev"],
            "noguchi_ratio_mean": noguchi_ratio["mean"],
            "noguchi_ratio_std": noguchi_ratio["std_dev"],
            "classification": noguchi_classification
    }


    return literature_by_id_noguchi

def compare_noguchi_hirata() -> dict:
    hirata = return_lit_id_hirata()
    noguchi = return_lit_id_noguchi()
    result = {}
    for id, crater_info in noguchi.items():
        nog_std = crater_info["noguchi_diameter_std"]
        nog_diam_mean = crater_info["noguchi_diameter_mean"]

        nog_min_diam = nog_diam_mean - nog_std
        nog_max_diam = nog_diam_mean + nog_std

        hirata_diam = hirata[id]["hirata_diameter_m"]
        mean_err = np.abs(nog_diam_mean - hirata_diam)
        gap_error = 0.0

        if hirata_diam < nog_min_diam:
            gap_error = nog_min_diam - hirata_diam
        elif hirata_diam > nog_max_diam:
            gap_error = hirata_diam - nog_max_diam


        result[id] = {
            "mean_err" : mean_err,
            "gap_error": gap_error,
            "hirata_diam": hirata_diam,
            "noguchi_diam": nog_diam_mean,
            "noguchi_std": nog_std
        }



    return result
