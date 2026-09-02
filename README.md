# APBT v2 crater-morphometry release

This repository is the reproducible computational release for the APBT v2 and
CIRCLE crater-morphometry comparison on Bennu, Ryugu, Itokawa, and Didymos. The
active scientific release is `final-all-manual-224`.

The repository contains the complete inputs, measurements, analysis code, and
derived outputs needed to validate the release. Manuscript source and working
documents are intentionally excluded.

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
three stages with one command:

```bash
python run_pipeline.py all --run complete --body all --method both
```

For a fast review using the completed measurement JSON files, add
`--use-saved` to `measure` or `all`:

```bash
python run_pipeline.py all --run quick-review --body Ryugu --method apbt \
  --reference hirata noguchi --use-saved
```

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
  manifest.json                 release counts and policies
  <Body>/
    seeds.json                  final active seed groups
    *_config.json               measurement parameters
    literature.json             literature comparison data
    meshes/*.vtk                meshes referenced by active seeds
results/
  release_manifest.json         release labels and integrity hashes
  <Body>/                       completed measurements and run records
  analysis/                     reproducibly generated tables and figures
methods/                        crater measurement algorithms
pipeline/                       output generation and literature comparisons
scripts/                        measurement entry points and path mapping
data/                           literature-data adapters
utils/                          shared geometry and identifier helpers
run_pipeline.py                 simple public pipeline entry point
```

The eight `results/<Body>/{apbt,circle}.json` files are the completed
measurements and the source of truth for derived analysis. Their integrity,
together with final inputs and the measurement implementation, is recorded in
`results/release_manifest.json`.
