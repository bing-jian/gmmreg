# GMM Registration — Experiment Scripts

These scripts reproduce the point cloud registration experiments from the gmmreg PAMI paper using the `gmmreg_demo` binary.

## Running with Docker (recommended)

A multi-stage `Dockerfile` at the repo root builds the C++ binary and sets up a Python 3.12 environment in one step.

```bash
# From the repo root — data paths default to ./data/dragon_stand and ./data/lounge
make run-dragon
make run-lounge

# Override data paths
make run-dragon DRAGON_DATA=/path/to/dragon_stand
make run-lounge LOUNGE_DATA=/path/to/lounge

# Pass extra script arguments
make run-dragon ARGS="--intervals 1,2 --voxel_size 0.003"

# Build image only
make build

# Drop into a shell inside the container for debugging
make shell

# Remove the image
make clean
```

> **Note:** `--visualize` flags require a display and will not work inside the container. Run visualizations locally instead.

---

## Prerequisites (running locally)

### 1. Build the binary

The scripts expect the compiled binary at `../C++/build/gmmreg_demo` (Linux/macOS) or `../C++/build/gmmreg_demo.exe` (Windows). Follow the build instructions in the top-level README before proceeding.

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

Python 3.12 is recommended. open3d does not support Python 3.13+.

---

## Experiment 1 — Stanford Dragon Stand

Reproduces Section 6.1 of the gmmreg PAMI paper. The dataset contains 15 scans of a dragon figurine taken at 24-degree intervals (0°, 24°, ..., 336°). The script runs all pairwise rigid registrations for pairs separated by 1, 2, and 3 steps (24°, 48°, 72°) and reports rotation errors.

### Get the dataset

Download the Stanford Dragon Stand dataset:
- URL: http://graphics.stanford.edu/data/3Dscanrep/
- Extract so that the PLY files and `dragonStandRight.conf` are under `../data/dragon_stand/`.

Expected layout:
```
../data/dragon_stand/
    dragonStandRight.conf
    dragonStandRight_0.ply
    dragonStandRight_24.ply
    ...
    dragonStandRight_336.ply
```

### Run

```bash
python dragon_expts.py
```

With all options explicit:

```bash
python dragon_expts.py \
    --data_dir ../data/dragon_stand \
    --config_tmpl ./dragon_stand.ini \
    --working_dir ./dragon_output \
    --voxel_size 0.005 \
    --intervals 1,2,3
```

To also visualize an example result after registration:

```bash
python dragon_expts.py --intervals 1 --visualize
```

### Key arguments

| Argument | Default | Description |
|---|---|---|
| `--data_dir` | `../data/dragon_stand` | Directory with PLY files and ground truth conf |
| `--config_tmpl` | `./dragon_stand.ini` | gmmreg config template |
| `--working_dir` | auto (UUID under `/tmp`) | Output directory for results |
| `--voxel_size` | `0.005` | Voxel size for downsampling PLY files |
| `--intervals` | `1,2,3` | Step sizes to test (1=24°, 2=48°, 3=72°, 4=96°) |
| `--gmmreg_exe` | auto-detected | Path to `gmmreg_demo` binary |
| `--visualize` | off | Show an example registration result with Open3D |

### Expected output

For each step size the script prints rotation error statistics across all 15 pairs:

```
Summary stats for registering pairs that are 24 degrees away from each other:
<avg_score, min_score, max_score, median_score>: ...
<# of high accuracy scores (>.99)>: N out of 15)
<avg_err, min_err, max_err, median_err>: ..., ..., ..., ... (in degrees)
<# of small errors (<3 degree)>: N out of 15)
```

---

## Experiment 2 — Stanford Lounge (RGB-D)

Tests pairwise rigid registration on consecutive RGB-D frames from an indoor scene. Depth images are converted to point clouds and registered with gmmreg.

### Get the dataset

```bash
wget https://github.com/isl-org/open3d_downloads/releases/download/20220301-data/LoungeRGBDImages.zip
unzip LoungeRGBDImages.zip -d ../data/lounge
```

Expected layout after extraction:
```
../data/lounge/
    lounge_trajectory.log
    depth/
        000001.png
        000002.png
        ...
        003000.png
```

### Run

To register and visualize a single pair:

```bash
python lounge_expts.py                              # uses default ../data/lounge
python lounge_expts.py --data_path /your/data/path
```

By default (see `__main__` block) this registers frames 1 and 11, then runs the full batch over all 2995 consecutive-5 pairs.

To change the pair or disable the full batch, edit the `__main__` block at the bottom of `lounge_expts.py`:

```python
# Register one pair with visualization
run_pairwise_registration(1, 11, visualize=True)

# Run full batch (2995 pairs, takes a while)
main()
```

### Expected output

```
Input pairs with pose difference ~= X.XX degrees
<avg_err, min_err, max_err, median_err>: ..., ..., ..., ... (in degrees)
<# of small errors (<3 degree)>: N out of 2995)
Registration time: <avg_time, min_time, max_time, median_time>: ... (in milliseconds)
```

Full batch results are saved to `./tmp/loung_expts_results.txt`.

---

## Output files

Both scripts write intermediate and final results under a `./tmp/` directory (created automatically):

| File | Contents |
|---|---|
| `tmp/model.txt` | Downsampled model point cloud (input to binary) |
| `tmp/scene.txt` | Downsampled scene point cloud (input to binary) |
| `tmp/final_rigid.txt` | Registration parameters `[q0 q1 q2 q3 tx ty tz]` |
| `tmp/final_rigid_matrix.txt` | 4×4 rigid transformation matrix |
| `tmp/transformed_model.txt` | Model after applying the estimated transform |
| `tmp/elapsed_time_in_ms.txt` | Core registration time reported by the binary |

For the dragon batch experiments, per-pair results are written to subdirectories under `--working_dir` named `{i}_vs_{j}/`.
