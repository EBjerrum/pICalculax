from setuptools import setup, find_packages

setup(
    name="picalculax",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "rdkit",
        "numpy",
    ],
    author="Esben Jannik Bjerrum",
    author_email="esben@cheminformania.com",  # Add if available
    description="A package for calculating pI values of peptides and modified peptides",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/EBjerrum/pICalculax",  # Update with actual repo URL
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "picalculax=picalculax.cli:main",
        ],
    },
) 