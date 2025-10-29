from setuptools import setup, find_packages

setup(
    name="cigar_bia",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pysam",
        "matplotlib"
    ],
    author="Your Name",
    description="CIGAR-Based Indels Analysis for KO/KI evaluation",
    url="https://github.com/yourusername/cigar_bia",
    python_requires=">=3.8"
)
