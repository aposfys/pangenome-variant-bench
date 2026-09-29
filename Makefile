.PHONY: install data run test lint clean all

PYTHON ?= python3
SAMPLE ?= HG002
REGION ?= chr20

all: test

install:
	$(PYTHON) -m pip install -e ".[dev]"

## Fetch the GIAB truth set, high-confidence BED and stratification BEDs for $(REGION)
data:
	$(PYTHON) -m panbench.fetch --sample $(SAMPLE) --region $(REGION)

## The planned caller comparison. main.nf is a skeleton that runs no process yet.
run:
	nextflow run . -profile docker,laptop --sample $(SAMPLE) --region $(REGION)

## The evaluation library only, with no Nextflow, containers or reference genome
test:
	$(PYTHON) -m pytest -q

lint:
	ruff check src tests && ruff format --check src tests && mypy src

clean:
	rm -rf work/ .nextflow* results/nextflow_*
	find . -name __pycache__ -type d -exec rm -rf {} +
