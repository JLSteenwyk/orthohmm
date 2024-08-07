profile:
	python3 -m cProfile -s 'time' orthohmm-runner.py ./tests/samples/ > large_profile.txt

run:
	python3 -m orthohmm-runner ./tests/samples/

run.simple:
	python3 -m orthohmm-runner ./tests/samples/	

install:
	# install so orthohmm command is available in terminal
	python3 setup.py install

develop:
	# https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode
	python3 setup.py develop

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
