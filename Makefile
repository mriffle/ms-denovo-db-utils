# Local development harness.
#
# CI runs these same targets (see .github/workflows/ci.yml) so that local and
# CI results cannot drift apart.

.DEFAULT_GOAL := help

PY      := python3
VENV    := .venv
BIN     := $(VENV)/bin
STAMP   := $(VENV)/.install-stamp
IMAGE   := ms-denovo-db-utils:dev

.PHONY: help install lint format typecheck test test-docker image check clean

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

$(STAMP): pyproject.toml
	$(PY) -m venv $(VENV)
	$(BIN)/python -m pip install --upgrade pip --quiet
	$(BIN)/python -m pip install -e '.[dev]' --quiet
	@touch $(STAMP)

install: $(STAMP) ## Create .venv and install the package plus dev tools

lint: install ## ruff check + format check + Dockerfile lint
	$(BIN)/ruff check .
	$(BIN)/ruff format --check .
	@if command -v hadolint >/dev/null 2>&1; then \
		hadolint Dockerfile; \
	else \
		echo "note: hadolint not installed locally; Dockerfile lint runs in CI"; \
	fi

format: install ## Auto-fix lint issues and format
	$(BIN)/ruff check --fix .
	$(BIN)/ruff format .

typecheck: install ## mypy --strict
	$(BIN)/mypy

test: install ## Unit and golden tests (no Docker)
	$(BIN)/pytest -m "not docker" --cov --cov-report=term-missing

image: ## Build the container image locally
	docker build -t $(IMAGE) .

test-docker: install ## Container smoke tests; skips cleanly if Docker is absent
	@if docker info >/dev/null 2>&1; then \
		$(MAKE) image && $(BIN)/pytest -m docker; \
	else \
		echo "SKIP: no usable Docker daemon; container tests not run"; \
	fi

check: lint typecheck test test-docker ## Everything CI runs

clean: ## Remove build artefacts and caches
	rm -rf $(VENV) .pytest_cache .mypy_cache .ruff_cache .coverage htmlcov \
	       build dist src/*.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
