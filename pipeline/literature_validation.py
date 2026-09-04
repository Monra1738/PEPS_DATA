"""Publication literature-comparison definitions.

These are analysis mappings, not gates that can reject or change measurements.
"""

from __future__ import annotations

from data.literature import return_lit_id_hirata, return_lit_id_noguchi
from data.literature_bennu import (
    return_lit_id_bierhaus,
    return_lit_id_daly,
    return_lit_id_deshapriya,
)
from data.literature_didymos import return_lit_id_barnouin
from data.literature_itokawa import return_lit_id_naru_hirata


CASES: list[dict] = []


def _add(
    body: str,
    method: str,
    metric: str,
    reference: str,
    literature_fn,
    value_key: str,
    std_key: str | None,
) -> None:
    CASES.append(
        {
            "asteroid": body,
            "method": method,
            "metric": metric,
            "reference": reference,
            "literature_fn": literature_fn,
            "literature_value_key": value_key,
            "literature_std_key": std_key,
        }
    )


for _method in ("APBT", "CIRCLE"):
    _add("Bennu", _method, "depth", "Daly", return_lit_id_daly,
         "daly_depth_mean", "daly_depth_std")
    _add("Bennu", _method, "diameter", "Daly", return_lit_id_daly,
         "daly_diameter_mean", "daly_diameter_std")
    _add("Bennu", _method, "diameter", "Bierhaus", return_lit_id_bierhaus,
         "bierhaus_diameter_mean", None)
    _add("Bennu", _method, "diameter", "Deshapriya", return_lit_id_deshapriya,
         "deshapriya_diameter_m", "deshapriya_diameter_std_m")
    _add("Ryugu", _method, "depth", "Noguchi", return_lit_id_noguchi,
         "noguchi_depth_mean", "noguchi_depth_std")
    _add("Ryugu", _method, "diameter", "Noguchi", return_lit_id_noguchi,
         "noguchi_diameter_mean", "noguchi_diameter_std")
    _add("Ryugu", _method, "diameter", "Hirata", return_lit_id_hirata,
         "hirata_diameter_m", None)
    _add("Itokawa", _method, "depth", "Naru-Hirata", return_lit_id_naru_hirata,
         "naru_hirata_depth_m", None)
    _add("Itokawa", _method, "diameter", "Naru-Hirata", return_lit_id_naru_hirata,
         "naru_hirata_diameter_m", "naru_hirata_diameter_std_m")
    _add("Didymos", _method, "diameter", "Barnouin", return_lit_id_barnouin,
         "Barnouin_diameter_m", "Barnouin_diameter_std_m")
