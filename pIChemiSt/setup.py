from setuptools import find_packages, setup
setup(
    name="pichemist",
    packages=find_packages(),
    version="0.0.1",
    author="Gian Marco Ghiandoni",
    author_email="ghiandoni.g@gmail.com",
    license="Apache",
    install_requires=["matplotlib", "rdkit"],
    setup_requires=["pytest-runner", "flake8"],
    tests_require=["pytest==7.1.2"],
    test_suite="test",
    package_dir={"pichemist": "pichemist"},
    url="https://github.com/AstraZeneca/peptide-tools/blob/"
        "master/pIChemiSt"
)
