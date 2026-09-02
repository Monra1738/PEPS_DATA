"""Body-specific paths shared by the APBT and CIRCLE batch entry points."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class BatchPaths:
    body: str
    mesh_dir: Path
    seeds: Path
    config: Path
    output: Path


APBT_PATHS = {
    "Bennu": BatchPaths(
        "Bennu",
        Path("inputs/Bennu/meshes"),
        Path("inputs/Bennu/seeds.json"),
        Path("inputs/Bennu/apbt_config.json"),
        Path("results/Bennu/apbt.json"),
    ),
    "Ryugu": BatchPaths(
        "Ryugu",
        Path("inputs/Ryugu/meshes"),
        Path("inputs/Ryugu/seeds.json"),
        Path("inputs/Ryugu/apbt_config.json"),
        Path("results/Ryugu/apbt.json"),
    ),
    "Itokawa": BatchPaths(
        "Itokawa",
        Path("inputs/Itokawa/meshes"),
        Path("inputs/Itokawa/seeds.json"),
        Path("inputs/Itokawa/apbt_config.json"),
        Path("results/Itokawa/apbt.json"),
    ),
    "Didymos": BatchPaths(
        "Didymos",
        Path("inputs/Didymos/meshes"),
        Path("inputs/Didymos/seeds.json"),
        Path("inputs/Didymos/apbt_config.json"),
        Path("results/Didymos/apbt.json"),
    ),
}


CIRCLE_PATHS = {
    "Bennu": BatchPaths(
        "Bennu",
        Path("inputs/Bennu/meshes"),
        Path("inputs/Bennu/seeds.json"),
        Path("inputs/Bennu/circle_config.json"),
        Path("results/Bennu/circle.json"),
    ),
    "Ryugu": BatchPaths(
        "Ryugu",
        Path("inputs/Ryugu/meshes"),
        Path("inputs/Ryugu/seeds.json"),
        Path("inputs/Ryugu/circle_config.json"),
        Path("results/Ryugu/circle.json"),
    ),
    "Itokawa": BatchPaths(
        "Itokawa",
        Path("inputs/Itokawa/meshes"),
        Path("inputs/Itokawa/seeds.json"),
        Path("inputs/Itokawa/circle_config.json"),
        Path("results/Itokawa/circle.json"),
    ),
    "Didymos": BatchPaths(
        "Didymos",
        Path("inputs/Didymos/meshes"),
        Path("inputs/Didymos/seeds.json"),
        Path("inputs/Didymos/circle_config.json"),
        Path("results/Didymos/circle.json"),
    ),
}


def project_path(path_like: str | Path) -> Path:
    """Return an absolute path, resolving relative paths from the project root."""

    path = Path(path_like)
    return path if path.is_absolute() else ROOT / path


def portable_path(path_like: str | Path) -> str:
    """Serialize project paths relative to the repository when possible."""

    path = Path(path_like).resolve()
    try:
        return path.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)
