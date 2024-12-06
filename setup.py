from setuptools import setup, find_packages

setup(
    name="anndata-validator",
    version="0.1.0",
    description="A package to validate AnnData objects.",
    author="Penny Cuda, HIVE CMU TC",
    license="GNU",
    packages=find_packages(),
    install_requires=[
        "anndata",
    ],
    python_requires=">=3.8",
)
