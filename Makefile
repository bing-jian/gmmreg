IMAGE_BUILDER ?= gmmreg-builder
IMAGE_TESTER  ?= gmmreg-tester
DATA_DIR      ?= $(PWD)/data

.PHONY: build-builder compile build-tester build run-dragon run-lounge shell clean

# Build the C++ build environment image (installs VXL and toolchain).
build-builder:
	docker build -f Dockerfile.builder -t $(IMAGE_BUILDER) .

# Compile the C++ binary by mounting the source tree into the builder image.
# Output lands in C++/build/ on the host.
compile: build-builder
	docker run --rm -v $(PWD)/C++:/workspace/C++ $(IMAGE_BUILDER)

# Build the Python test environment image (depends on builder image).
build-tester: build-builder
	docker build -f Dockerfile.tester -t $(IMAGE_TESTER) .

# Build both images (does not compile C++; run 'make compile' separately).
build: build-tester

run-dragon: build-tester
	docker run --rm \
		-v $(PWD)/expts:/workspace/expts \
		-v $(PWD)/C++/build:/workspace/C++/build:ro \
		-v $(DATA_DIR):/workspace/data:ro \
		$(IMAGE_TESTER) \
		python dragon_expts.py --data_dir /workspace/data/dragon_stand $(ARGS)

run-lounge: build-tester
	docker run --rm \
		-v $(PWD)/expts:/workspace/expts \
		-v $(PWD)/C++/build:/workspace/C++/build:ro \
		-v $(DATA_DIR):/workspace/data:ro \
		$(IMAGE_TESTER) \
		python lounge_expts.py --data_path /workspace/data/lounge $(ARGS)

shell: build-tester
	docker run --rm -it \
		-v $(PWD)/expts:/workspace/expts \
		-v $(PWD)/C++/build:/workspace/C++/build \
		-v $(DATA_DIR):/workspace/data \
		$(IMAGE_TESTER) bash

clean:
	docker rmi $(IMAGE_BUILDER) $(IMAGE_TESTER)
