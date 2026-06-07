# Docker Build and Test Workflow

Two Dockerfiles cleanly separate the C++ build environment from the Python test environment.

## Images

| Dockerfile | Image | Purpose |
|---|---|---|
| `Dockerfile.builder` | `gmmreg-builder` | Ubuntu 24.04 + CMake + VXL; compiles the C++ binary |
| `Dockerfile.tester` | `gmmreg-tester` | Extends builder with Python 3.12 and experiment dependencies |

The `C++/` source tree and `data/` are never baked into images — they are always mounted at runtime. This means you can edit source or swap datasets without rebuilding images.

---

## Quickstart

```bash
# 1. Build both images (one-time, or after changing Dockerfiles)
make build

# 2. Compile the C++ binary (re-run whenever C++ source changes)
make compile

# 3. Run experiments
make run-dragon
make run-lounge
```

---

## Make targets

| Target | Description |
|---|---|
| `make build-builder` | Build `gmmreg-builder` image (VXL + toolchain + GoogleTest) |
| `make compile` | Mount `C++/` into the builder and produce the binary in `C++/build/` |
| `make test` | Build and run all GoogleTest/CTest unit tests inside the builder container |
| `make shell-builder` | Interactive bash shell in the builder container with `C++/`, `expts/`, `data/` mounted |
| `make build-tester` | Build `gmmreg-tester` image (depends on `gmmreg-builder`) |
| `make build` | Shortcut: `build-builder` + `build-tester` |
| `make run-dragon` | Run the Dragon Stand experiment inside the tester container |
| `make run-lounge` | Run the Lounge RGB-D experiment inside the tester container |
| `make shell` | Open an interactive bash session in the tester container |
| `make clean` | Remove both images |

### Overridable variables

| Variable | Default | Description |
|---|---|---|
| `IMAGE_BUILDER` | `gmmreg-builder` | Builder image tag |
| `IMAGE_TESTER` | `gmmreg-tester` | Tester image tag |
| `DATA_DIR` | `$(PWD)/data` | Host path mounted as `/workspace/data` |
| `ARGS` | _(empty)_ | Extra arguments forwarded to the Python script |

Examples:

```bash
# Use a different data directory
make run-dragon DATA_DIR=/mnt/datasets

# Pass extra flags to the script
make run-dragon ARGS="--intervals 1,2 --voxel_size 0.003"

# Override image names (e.g. for a registry push)
make build IMAGE_BUILDER=myrepo/gmmreg-builder IMAGE_TESTER=myrepo/gmmreg-tester
```

---

## Volume mounts at runtime

| Host path | Container path | Used by |
|---|---|---|
| `C++/` | `/workspace/C++` | `compile` (builder) |
| `C++/build/` | `/workspace/C++/build` | `run-*`, `shell` (tester) |
| `expts/` | `/workspace/expts` | `run-*`, `shell` (tester) |
| `data/` | `/workspace/data` | `run-*`, `shell` (tester) |

Mounting `expts/` means you can edit Python scripts on the host and re-run without rebuilding the tester image.

---

## Rebuilding after changes

| What changed | Command needed |
|---|---|
| `Dockerfile.builder` or VXL | `make build-builder` then `make compile` then `make build-tester` |
| C++ source | `make compile` |
| C++ tests | `make test` |
| `Dockerfile.tester` or `requirements.txt` | `make build-tester` |
| Python scripts in `expts/` | Nothing — changes are live via the mount |
