# ── Stage 1: Build C++ binary ─────────────────────────────────────────────
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Build VXL (provides vnl/vnl_algo/vcl)
RUN git clone --depth 1 https://github.com/vxl/vxl.git /vxl-src && \
    cmake -B /vxl-build -S /vxl-src \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=OFF \
        -DBUILD_TESTING=OFF \
        -DBUILD_DOCUMENTATION=OFF \
        -DVXL_BUILD_EXAMPLES=OFF && \
    cmake --build /vxl-build --parallel $(nproc)

ENV VXL_DIR=/vxl-build

COPY C++/ /workspace/C++/
WORKDIR /workspace/C++
RUN cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build --parallel $(nproc)


# ── Stage 2: Python test environment ──────────────────────────────────────
FROM python:3.12-slim AS runner

# System libraries required by open3d and opencv-python
RUN apt-get update && apt-get install -y \
    libgomp1 \
    libgl1 \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

# Copy only the compiled binary into the path the Python scripts expect
# (expts/ scripts resolve BINARY_DIR as '../C++/build' relative to /workspace/expts)
COPY --from=builder /workspace/C++/build/gmmreg_demo /workspace/C++/build/gmmreg_demo

# Install Python dependencies (expts/ is mounted at runtime, not copied)
COPY expts/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --root-user-action=ignore -r /tmp/requirements.txt

# Mount expts/ at runtime: docker run -v $(PWD)/expts:/workspace/expts
WORKDIR /workspace/expts
