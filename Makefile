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
	python3 -m pytest tests/unit

test.integration:
	python3 -m pytest tests/integration

test.fast:
	python3 -m pytest tests/unit -m "not slow"
	python3 -m pytest tests/integration -m "not slow"

# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python3 -m pytest tests/unit --cov=./ --cov-report=xml:unit.coverage.xml

coverage.integration:
	python3 -m pytest tests/integration --cov=./ --cov-report=xml:integration.coverage.xml
