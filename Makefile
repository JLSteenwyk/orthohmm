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
	python3 -m pytest -m "not integration"

test.integration:
	rm -rf ./tests/samples/orthohmm_*
	python3 -m pytest --basetemp=output -m "integration"

test.fast:
	python3 -m pytest -m "not (integration or slow)"
	rm -rf ./tests/samples/orthohmm_*
	python3 -m pytest --basetemp=output -m "integration and not slow"

# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python3 -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm -rf ./tests/samples/orthohmm_*
	python3 -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml
