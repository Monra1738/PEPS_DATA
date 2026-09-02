# APBT v2 crater-morphometry release

This repository is the reproducible computational release for the APBT v2 and
CIRCLE crater-morphometry comparison on Bennu, Ryugu, Itokawa, and Didymos. The
active scientific release is `final-all-manual-224`.

The repository contains the complete inputs and analysis code needed to run the
project from scratch. Generated results and manuscript source are intentionally
excluded.

## Release snapshot

| Body | Seeded craters | Final seed groups | Valid APBT | Valid CIRCLE |
|---|---:|---:|---:|---:|
| Bennu | 45 | 135 | 135 | 123 |
| Ryugu | 77 | 230 | 230 | 229 |
| Itokawa | 21 | 224 | 217 | 158 |
| Didymos | 20 | 60 | 60 | 53 |
| **Total** | **163** | **649** | **642** | **563** |

Mesh coordinates are supplied in kilometres. Saved diameters and depths are in
metres. The 163 VTK files under `inputs/` are required measurement inputs and
are stored as ASCII mesh data.

## Installation

Create the pinned Conda environment:

```bash
conda env create -f environment.yml
conda activate apbt-repro
```

## Run the project

The pipeline has three simple stages. Give the run a name so its files remain
separate under `runs/<name>/`.

```bash
# 1. Calculate one or more measurements
python run_pipeline.py measure --run ryugu-test --body Ryugu --method apbt

# 2. Compare those results with selected literature
python run_pipeline.py compare --run ryugu-test --reference hirata noguchi

# 3. Create the final CSV tables and PNG figures
python run_pipeline.py outputs --run ryugu-test
```

Use `--body all` and `--method both` for the complete project. To perform all
stages and build the full publication analysis with one command:

```bash
python run_pipeline.py all --run complete --body all --method both
```

The `all` command calculates measurements, performs literature comparisons,
creates run-level outputs under `runs/complete/outputs/`, and automatically
builds the publication tables and figures under `results/analysis/`. The
publication analysis is available only for all four bodies and both methods.

## Itokawa population

`inputs/Itokawa/seeds.json` contains all 224 structurally valid manual
four-index groups across 21 craters. The original manual pool contained 226
groups; two groups were excluded before measurement because each repeated a
vertex index. No random or replacement mesh points were generated.

APBT saved 217 finite positive results and rejected seven invalid rim/floor
geometries. CIRCLE saved 158 finite positive results and rejected 66. There are
156 seed identities with a valid result from both methods. The run-status JSON
files preserve every attempt.

## Repository layout

```text
environment.yml                 pinned Python environment
inputs/
  <Body>/
    seeds.json                  final active seed groups
    *_config.json               measurement parameters
    literature.json             literature comparison data
    meshes/*.vtk                meshes referenced by active seeds
methods/                        crater measurement algorithms
pipeline/                       output generation and literature comparisons
scripts/                        measurement entry points and path mapping
data/                           literature-data adapters
utils/                          shared geometry and identifier helpers
run_pipeline.py                 simple public pipeline entry point
```

The `results/` directory is created locally when you run measurements and
analysis. It is intentionally not included in this repository.
