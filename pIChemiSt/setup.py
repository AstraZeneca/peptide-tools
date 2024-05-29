from setuptools import find_packages
from setuptools import setup

__doc__ = (
    "Calculates the isoelectric point (pI) of a peptide/protein "
    "using the Henderson-Hasselbalch equation with the predicted "
    "pKa values of its amino acids."
    ""
)

# read the contents of your README file
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="pichemist",
    packages=find_packages(),
    version="0.1.0",
    description=__doc__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Gian Marco Ghiandoni",
    author_email="ghiandoni.g@gmail.com",
    license="Apache 2.0",
    install_requires=["matplotlib", "rdkit"],
    setup_requires=["pytest-runner", "flake8"],
    tests_require=["pytest==7.1.2"],
    test_suite="test",
    package_dir={"pichemist": "pichemist"},
    include_package_data=True,
    url="https://github.com/AstraZeneca/peptide-tools/blob/" "master/pIChemiSt",
    entry_points={
        "console_scripts": ["pichemist=pichemist.cli:main"],
    },
)
