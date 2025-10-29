from setuptools import setup, find_packages
from cigar_bia import __version__

setup(
    name="cigar_bia",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pysam",
        "matplotlib"
    ],
    author="Francesco Patanè",
    description="CIGAR-Based Indels Analysis for KO/KI evaluation",
    url="https://github.com/francescopatane96/cigar_bia",
    python_requires=">=3.8"
)
