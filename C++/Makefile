BUILDER_IMAGE ?= gmmreg-builder

.PHONY: builder shell

builder:
	docker buildx build -t $(BUILDER_IMAGE) -f Dockerfile.builder .

shell:
	docker run --rm -it -v $(PWD):/workspace $(BUILDER_IMAGE) bash
