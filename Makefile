PYTEST_ARGS ?=
TEST_RECEIPT_DIR ?= .
KERNEL_BUILD_DIR ?= .cache/checkout-kernels

.PHONY: install develop kernels test test.unit test.integration test.fast test.coverage coverage.unit coverage.integration coverage.combine

install:
	python3 -m pip install .

develop:
	python3 -m pip install -e .

kernels:
	python3 scripts/build_ci_checkout_kernels.py --root "$(CURDIR)" --output "$(abspath $(KERNEL_BUILD_DIR))"

test: test.unit test.integration

test.unit:
	python3 -m pytest tests --ignore=tests/integration --junitxml="$(TEST_RECEIPT_DIR)/unit.xml" $(PYTEST_ARGS)

test.integration:
	python3 -m pytest tests/integration --junitxml="$(TEST_RECEIPT_DIR)/integration.xml" $(PYTEST_ARGS)

test.fast:
	$(MAKE) test PYTEST_ARGS='-m "not slow" $(PYTEST_ARGS)'

test.coverage: coverage.unit coverage.integration
	$(MAKE) coverage.combine

coverage.unit:
	COVERAGE_FILE="$(abspath .coverage-unit)" python3 -m pytest tests --ignore=tests/integration --cov --cov-config=.coveragerc --cov-report=xml:unit.coverage.xml --junitxml="$(TEST_RECEIPT_DIR)/unit.xml" $(PYTEST_ARGS)

coverage.integration:
	COVERAGE_FILE="$(abspath .coverage-integration)" python3 -m pytest tests/integration --cov --cov-config=.coveragerc --cov-report=xml:integration.coverage.xml --junitxml="$(TEST_RECEIPT_DIR)/integration.xml" $(PYTEST_ARGS)

# Stage files are not .coverage.* parallel fragments and remain auditable.
coverage.combine:
	python3 -m coverage combine --keep .coverage-unit .coverage-integration
	python3 -m coverage xml -o coverage.xml
	python3 -m coverage json -o coverage.json
	python3 -m coverage report --include="orthohmm/*" --fail-under=85
