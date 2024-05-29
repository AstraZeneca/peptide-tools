.PHONY: build install test coverage

build:
	python setup.py bdist_wheel
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

upload_test:
	twine upload -r testpypi dist/*

upload:
	twine upload -r pypi dist/*

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
