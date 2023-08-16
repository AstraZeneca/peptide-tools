.PHONY: build install test coverage

build:
	python setup.py sdist

install:
	pip install dist/pichemist-*

test_core:
	pytest -vvs test

test_acd:
	pytest -vvs test -m acd

test:
	make test_core
	make test_acd

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
