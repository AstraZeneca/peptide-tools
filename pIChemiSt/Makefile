.PHONY: build install test coverage

build:
	python setup.py bdist_wheel

install:
	pip install dist/cse-*.whl

test_core:
	pytest -vvs test

test_acd:
	pytest -vvs test -m acd

test_full:
	make test_core
	make test_acd

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
