.PHONY: build install test coverage

build:
	python setup.py bdist_wheel

install:
	pip install dist/cse-*.whl

test:
	pytest -vvs test
	# pytest -vs test/test_api.py

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
