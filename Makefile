IMAGE   ?= gmmreg-expts
DRAGON_DATA ?= $(PWD)/data/dragon_stand
LOUNGE_DATA ?= $(PWD)/data/lounge

.PHONY: build run-dragon run-lounge shell clean

build:
	docker build -t $(IMAGE) .

run-dragon: build
	docker run --rm \
		-v $(PWD)/expts:/workspace/expts \
		-v $(DRAGON_DATA):/workspace/data/dragon_stand:ro \
		$(IMAGE) \
		python dragon_expts.py --data_dir /workspace/data/dragon_stand $(ARGS)

run-lounge: build
	docker run --rm \
		-v $(PWD)/expts:/workspace/expts \
		-v $(LOUNGE_DATA):/workspace/data/lounge:ro \
		$(IMAGE) \
		python lounge_expts.py --data_path /workspace/data/lounge $(ARGS)

shell: build
	docker run --rm -it \
		-v $(PWD)/expts:/workspace/expts \
		-v $(DRAGON_DATA):/workspace/data/dragon_stand:ro \
		-v $(LOUNGE_DATA):/workspace/data/lounge:ro \
		$(IMAGE) bash

clean:
	docker rmi $(IMAGE)
