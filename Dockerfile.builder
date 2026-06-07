FROM ubuntu:24.04

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

# C++ source is mounted at runtime, not copied.
# Run: docker run --rm -v $(PWD)/C++:/workspace/C++ gmmreg-builder
WORKDIR /workspace/C++
# Remove stale CMakeCache.txt and CMakeFiles/ to avoid path mismatch when the
# source tree was previously configured outside this container.
CMD rm -rf build/CMakeCache.txt build/CMakeFiles && \
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build --parallel $(nproc)
