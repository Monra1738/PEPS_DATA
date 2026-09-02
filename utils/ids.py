from pathlib import Path
import re


def extract_crater_id(name : str) -> int:
    stem = Path(name).stem

    pattern = r'\d+'
    m = re.search(pattern, stem)

    id = int(m.group(0))

    return id
