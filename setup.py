from setuptools import find_packages, setup

setup(
    name="checkv",
    version="0.1.0",
    packages=find_packages(),
    license="Modified BSD",
    description="Assess the quality of metagenome-assembled viral genomes.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "biopython",
        "importlib-metadata>=0.12; python_version<'3.8'",
        "numpy",
        "psutil",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["checkv=checkv.cli:cli"]},
    url="https://bitbucket.org/berkeleylab/checkv",
    keywords=[
        "bioinformatics",
        "genomics",
        "metagenomics",
        "viromics"
    ],
    author="Stephen Nayfach, Antonio Pedro Camargo, Simon Roux",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
