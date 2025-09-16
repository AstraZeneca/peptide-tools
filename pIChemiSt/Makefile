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
	version=$$(grep -Po "__version__\s*=\s*['\"]\K[^'\"]+" setup.py); \
	twine upload -r testpypi dist/pichemist-$$version* --verbose

upload:
	version=$$(grep -Po "__version__\s*=\s*['\"]\K[^'\"]+" setup.py); \
	twine upload -r pypi dist/pichemist-$$version* --verbose

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing

create_git_tag:
	version=$$(grep -Po "__version__\s*=\s*['\"]\K[^'\"]+" setup.py); \
	git tag v$$version
	git push --tags
