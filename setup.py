from setuptools import find_packages, setup

setup(
    name="checkv",
    version="0.1.0",
    packages=find_packages(),
    license="GNU General Public License v3.0",
    description="Identify and remove incorrectly binned contigs from metagenome-assembled genomes.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "biopython",
        "numpy",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["checkv=checkv.cli:cli"]},
    url="https://bitbucket.org/berkeleylab/checkv",
    keywords=[
        "bioinformatics",
        "genomics",
        "metagenomics"
    ],
    author="Stephen Nayfach, Antonio Pedro Camargo",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python :: 3",
    ],
)