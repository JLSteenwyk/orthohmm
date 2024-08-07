from os import path
from setuptools import setup

from orthohmm.version import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering",
]

REQUIRES = ["biopython>=1.84", "numpy>=1.24.0", "pandas>=2.2.2", "cython"]

setup(
    name="orthohmm",
    description="HMM-based inference of orthologous groups.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/orthohmm",
    packages=["orthohmm"],
    classifiers=CLASSIFIERS,
    entry_points={"console_scripts": ["orthohmm = orthohmm.orthohmm:main"]},
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
)

## push new version to pypi
# rm -rf dist
# python3 setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
