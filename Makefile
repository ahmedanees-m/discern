.PHONY: help install lint test data fmt

help:
	@echo "install  - create the conda env (environment.yml)"
	@echo "lint     - ruff check"
	@echo "fmt      - black + ruff --fix"
	@echo "test     - ruff check + pytest"
	@echo "data     - download/verify reference data (writes data/manifest.json)"

install:
	conda env create -f environment.yml

lint:
	ruff check .

fmt:
	black .
	ruff check --fix .

test:
	ruff check .
	pytest -q

data:
	python -m data.sources.build_manifest
