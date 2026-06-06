FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Build and install VXL (provides vnl/vnl_algo/vcl)
RUN git clone --depth 1 https://github.com/vxl/vxl.git /vxl-src && \
    cmake -B /vxl-build -S /vxl-src \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=OFF \
        -DBUILD_TESTING=OFF \
        -DBUILD_DOCUMENTATION=OFF \
        -DVXL_BUILD_EXAMPLES=OFF && \
    cmake --build /vxl-build --parallel $(nproc)

ENV VXL_DIR=/vxl-build

WORKDIR /workspace
