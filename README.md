# APBT v2 crater-morphometry release

This repository is the reproducible computational release for the APBT v2 and
CIRCLE crater-morphometry comparison on Bennu, Ryugu, Itokawa, and Didymos.
The current working release uses the cleaned Itokawa seed set described below.

It contains the measurement inputs, analysis code, and regenerated measurement
results used by the current analysis. Manuscript sources and other local
working documents are intentionally excluded.

## Current results

| Body | Craters | Active seed groups | Valid APBT | Valid CIRCLE |
|---|---:|---:|---:|---:|
| Bennu | 45 | 135 | 135 | 124 |
| Ryugu | 77 | 230 | 230 | 229 |
| Itokawa | 21 | 202 | 202 | 145 |
| Didymos | 20 | 60 | 60 | 54 |
| **Total** | **163** | **627** | **627** | **552** |

APBT completed successfully for all 627 active seed groups. CIRCLE produced
552 valid results; its 75 failed attempts are retained in the run-status JSON
files. Mesh coordinates are supplied in kilometres. Saved diameters and depths
are in metres. The 163 VTK files under `inputs/` are required measurement
inputs and are stored as ASCII mesh data.

## Installation

Create the pinned Conda environment:

```bash
conda env create -f environment.yml
conda activate apbt-repro
```

## Run the project

Run the complete analysis from the current inputs with a new run name:

```bash
python run_pipeline.py all \
  --run complete \
  --body all \
  --method both
```

The pipeline measures the selected bodies and methods, compares the results
with literature, and creates CSV summaries and PNG figures. Run-specific files
are written under `runs/<name>/outputs/`. Use a new run name when rerunning;
existing measurement files are kept by the workflow.

Individual stages can also be run separately:

```bash
# Calculate one or more measurements
python run_pipeline.py measure --run ryugu-test --body Ryugu --method apbt

# Compare results with selected literature references
python run_pipeline.py compare --run ryugu-test --reference hirata noguchi

# Create CSV tables and PNG comparison figures
python run_pipeline.py outputs --run ryugu-test
```

`--body all` selects all four bodies and `--method both` selects APBT and
CIRCLE. Valid run names contain only letters, numbers, `.`, `-`, and `_`.

## Verified analysis results

The regenerated report at `results/analysis/analysis_report.md` contains 20
standalone literature-validation cases:

| Body | Method | Metric | Reference | N | MAE (m) | MAPE (%) |
|---|---|---|---|---:|---:|---:|
| Bennu | APBT | depth | Daly | 33 | 1.340 | 28.886 |
| Bennu | APBT | diameter | Daly | 33 | 6.590 | 10.842 |
| Bennu | APBT | diameter | Bierhaus | 44 | 6.186 | 13.210 |
| Bennu | APBT | diameter | Deshapriya | 45 | 6.686 | 13.332 |
| Ryugu | APBT | depth | Noguchi | 77 | 0.916 | 17.978 |
| Ryugu | APBT | diameter | Noguchi | 77 | 5.569 | 11.482 |
| Ryugu | APBT | diameter | Hirata | 77 | 7.568 | 22.306 |
| Itokawa | APBT | depth | Naru-Hirata | 21 | 2.690 | 66.606 |
| Itokawa | APBT | diameter | Naru-Hirata | 21 | 4.976 | 9.657 |
| Didymos | APBT | diameter | Barnouin | 16 | 11.682 | 12.247 |
| Bennu | CIRCLE | depth | Daly | 32 | 1.966 | 39.058 |
| Bennu | CIRCLE | diameter | Daly | 32 | 8.128 | 14.175 |
| Bennu | CIRCLE | diameter | Bierhaus | 41 | 6.903 | 18.193 |
| Bennu | CIRCLE | diameter | Deshapriya | 42 | 8.019 | 19.139 |
| Ryugu | CIRCLE | depth | Noguchi | 77 | 1.123 | 21.574 |
| Ryugu | CIRCLE | diameter | Noguchi | 77 | 6.746 | 14.226 |
| Ryugu | CIRCLE | diameter | Hirata | 77 | 9.294 | 30.813 |
| Itokawa | CIRCLE | depth | Naru-Hirata | 19 | 3.719 | 83.562 |
| Itokawa | CIRCLE | diameter | Naru-Hirata | 19 | 7.301 | 13.458 |
| Didymos | CIRCLE | diameter | Barnouin | 14 | 17.531 | 15.982 |

The validation data are comparisons only; they do not reject or alter the
declared measurements.

## Itokawa population

`inputs/Itokawa/seeds.json` contains 202 active manual four-index groups across
21 craters. The original manual pool contained 226 groups. Two groups were
excluded because each repeated a vertex index, and 22 additional groups were
removed during APBT validation. No random or replacement mesh points were
generated.

APBT produced 202/202 valid results with no remaining invalid attempts. CIRCLE
produced 145/202 valid results; 57 attempts failed because of non-positive
rim/floor depth estimates. The run-status JSON files preserve every attempt.

## Repository layout

```text
environment.yml                 pinned Python environment
inputs/
  <Body>/
    seeds.json                  active seed groups
    *_config.json               measurement parameters
    literature.json             literature comparison data
    meshes/*.vtk                meshes referenced by active seeds
methods/                        crater measurement algorithms
pipeline/                       workflow, summaries, and analysis
scripts/                        measurement entry points and path mapping
data/                           literature-data adapters
run_pipeline.py                 public pipeline entry point
results/<Body>/                 regenerated measurement JSON files
results/analysis/               verified tables, figures, and report
runs/<name>/                    local staged-run outputs (created at runtime)
```
