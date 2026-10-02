PYTEST_ARGS ?=
TEST_RECEIPT_DIR ?= .

profile:
	python3 -m cProfile -s 'time' orthohmm-runner.py ./tests/samples/three_proteome_set/ |tee large_profile.txt

run:
	python3 -m orthohmm-runner ./tests/samples/

run.simple:
	python3 -m orthohmm-runner ./tests/samples/	

install:
	# install so orthohmm command is available in terminal.
	# pip (not `setup.py install`) so transitive deps like leidenalg /
	# python-igraph resolve to prebuilt wheels instead of compiling
	# from source — easy_install would otherwise fall back to building
	# from sdist and fail on hosts without igraph C headers.
	python3 -m pip install .

develop:
	python3 -m pip install -e .

test: test.unit test.integration

test.unit:
	python3 -m pytest tests --ignore=tests/integration $(PYTEST_ARGS)

test.integration:
	python3 -m pytest tests/integration $(PYTEST_ARGS)

test.fast:
	python3 -m pytest tests --ignore=tests/integration -m "not slow" $(PYTEST_ARGS)
	python3 -m pytest tests/integration -m "not slow" $(PYTEST_ARGS)

# Complete coverage profile, including separately supplied raw-source cases.
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python3 -m pytest tests --ignore=tests/integration --cov=./ --cov-report=xml:unit.coverage.xml $(PYTEST_ARGS)

coverage.integration:
	python3 -m pytest tests/integration --cov=./ --cov-report=xml:integration.coverage.xml $(PYTEST_ARGS)

# Public CI is not the complete raw-benchmark regression gate.
test.public.fast:
	python3 -m pytest tests --ignore=tests/integration -m "not slow and not raw_benchmark" --junitxml="$(TEST_RECEIPT_DIR)/public-unit.xml"
	python3 -m pytest tests/integration -m "not slow and not raw_benchmark" --junitxml="$(TEST_RECEIPT_DIR)/public-integration.xml"

test.public.coverage: coverage.public.unit coverage.public.integration

coverage.public.unit:
	python3 -m pytest tests --ignore=tests/integration -m "not raw_benchmark" --cov=./ --cov-report=xml:unit.coverage.xml --junitxml="$(TEST_RECEIPT_DIR)/public-unit.xml"

coverage.public.integration:
	python3 -m pytest tests/integration -m "not raw_benchmark" --cov=./ --cov-report=xml:integration.coverage.xml --junitxml="$(TEST_RECEIPT_DIR)/public-integration.xml"
