from pathlib import Path

from setuptools import find_packages, setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="hubmap-mudata-validator",
    version="1.0.7",
    description="A package to validate MuData objects for HuBMAP EPICs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Penny Cuda",
    author_email="pcuda@andrew.cmu.edu",
    license="GPLv3",
    packages=find_packages(),
    install_requires=[
        "muon",
    ],
    url="https://github.com/hubmapconsortium/mudata-validator",
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
