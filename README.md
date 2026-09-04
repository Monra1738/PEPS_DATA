# APBT v2 crater-morphometry release

This repository is the reproducible computational release for the APBT v2 and
CIRCLE crater-morphometry comparison on Bennu, Ryugu, Itokawa, and Didymos.
The active scientific release is `final-all-manual-224`.

It contains the measurement inputs, analysis code, and the declared final
measurement results used by the presented-paper analysis. Manuscript sources
and other local working documents are intentionally excluded.

## Release snapshot

| Body | Craters | Final seed groups | Valid APBT | Valid CIRCLE |
|---|---:|---:|---:|---:|
| Bennu | 45 | 135 | 135 | 123 |
| Ryugu | 77 | 230 | 230 | 229 |
| Itokawa | 21 | 224 | 217 | 158 |
| Didymos | 20 | 60 | 60 | 53 |
| **Total** | **163** | **649** | **642** | **563** |

Mesh coordinates are supplied in kilometres. Saved diameters and depths are
in metres. The 163 VTK files under `inputs/` are required measurement inputs
and are stored as ASCII mesh data.

## Installation

Create the pinned Conda environment:

```bash
conda env create -f environment.yml
conda activate apbt-repro
```

## Reproduce the release

The checked-in JSON files under `results/<Body>/` are the declared final
measurements. Use `--use-saved` to build a complete run from those files
without recalculating the measurements:

```bash
python run_pipeline.py all \
  --run complete \
  --body all \
  --method both \
  --use-saved
```

This performs measurement staging, literature comparison, and output
generation. It writes run-specific files to `runs/complete/outputs/` and
regenerates the publication tables and figures under `results/analysis/`.

To recalculate measurements from the inputs instead, omit `--use-saved`.
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

The checked-in report at `results/analysis/analysis_report.md` contains 18
standalone literature-validation cases. The aggregate error results are:

| Body | Method | Metric | Reference | N | MAE (m) | MAPE (%) |
|---|---|---|---|---:|---:|---:|
| Bennu | APBT | depth | Daly | 33 | 1.343 | 28.938 |
| Bennu | APBT | diameter | Daly | 33 | 6.504 | 10.757 |
| Bennu | APBT | diameter | Bierhaus | 44 | 6.251 | 13.266 |
| Bennu | APBT | diameter | Deshapriya | 45 | 6.749 | 13.387 |
| Ryugu | APBT | depth | Noguchi | 77 | 0.916 | 17.978 |
| Ryugu | APBT | diameter | Noguchi | 77 | 5.569 | 11.482 |
| Itokawa | APBT | depth | Naru-Hirata | 21 | 2.663 | 66.688 |
| Itokawa | APBT | diameter | Naru-Hirata | 21 | 5.104 | 9.832 |
| Didymos | APBT | diameter | Barnouin | 16 | 11.597 | 12.215 |
| Bennu | CIRCLE | depth | Daly | 32 | 1.948 | 38.803 |
| Bennu | CIRCLE | diameter | Daly | 32 | 8.003 | 14.051 |
| Bennu | CIRCLE | diameter | Bierhaus | 41 | 7.001 | 18.278 |
| Bennu | CIRCLE | diameter | Deshapriya | 42 | 8.115 | 19.222 |
| Ryugu | CIRCLE | depth | Noguchi | 77 | 1.123 | 21.574 |
| Ryugu | CIRCLE | diameter | Noguchi | 77 | 6.746 | 14.226 |
| Itokawa | CIRCLE | depth | Naru-Hirata | 20 | 4.143 | 82.063 |
| Itokawa | CIRCLE | diameter | Naru-Hirata | 20 | 7.527 | 13.257 |
| Didymos | CIRCLE | diameter | Barnouin | 14 | 17.854 | 16.102 |

The validation data are comparisons only; they do not reject or alter the
declared measurements.

## Itokawa population

`inputs/Itokawa/seeds.json` contains all 224 structurally valid manual
four-index groups across 21 craters. The original manual pool contained 226
groups; two groups were excluded before measurement because each repeated a
vertex index. No random or replacement mesh points were generated.

APBT saved 217 finite positive results and rejected seven invalid rim/floor
geometries. CIRCLE saved 158 finite positive results and rejected 66. There
are 156 seed identities with a valid result from both methods. The run-status
JSON files preserve every attempt.

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
pipeline/                       workflow, summaries, and analysis
scripts/                        measurement entry points and path mapping
data/                           literature-data adapters
run_pipeline.py                 public pipeline entry point
results/<Body>/                 declared final measurement JSON files
results/analysis/               verified tables, figures, and report
runs/<name>/                    local staged-run outputs (created at runtime)
```

```bash
for body in Bennu Ryugu Itokawa Didymos; do
    mkdir -p "results/$body"
    cp "runs/complete/measurements/$body/apbt.json" "results/$body/apbt.json"
    cp "runs/complete/measurements/$body/circle.json" "results/$body/circle.json"
  done

  python pipeline/analysis.py

  ```