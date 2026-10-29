.PHONY: help image test fixture run run-ieeg clean

SERVICE_NAME ?= processor-seizure-consensus
PYTHON ?= python3
IMAGE ?= $(SERVICE_NAME):dev
# A folder with ieeg/ (ieeg-prep's output) and methods/ (dynasd, ei, fragility outputs).
IN ?= data/real
export DOCKER_DEFAULT_PLATFORM = linux/amd64

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test     - pytest: synthetic fixture through the processor, both input layouts, Lambda handler (no Docker)"
	@echo "make fixture  - write the synthetic fixture to data/fixture/{linked,flat}/"
	@echo "make run      - build the image and run it in Docker over data/fixture/linked/, writing to data/output/"
	@echo "make run-ieeg - run the image over IN=$(IN) (ieeg/ + methods/), writing to data/output-real/"
	@echo "make clean    - remove data/output*/ and data/fixture/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

image:
	docker build --platform linux/amd64 -t $(IMAGE) .

test:
	$(PYTHON) -m pytest -q tests

fixture:
	$(PYTHON) -m tests.make_fixture data/fixture

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

run-ieeg: image
	rm -rf data/output-real && mkdir -p data/output-real
	docker run --rm --platform linux/amd64 --memory 1g --env-file dev.env \
		-v "$(abspath $(IN))":/data/input:ro -v "$(CURDIR)/data/output-real":/data/output $(IMAGE)

clean:
	rm -rf data/output data/output-real data/fixture
