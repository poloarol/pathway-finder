import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="pathway-finder-polo316",
    version="0.0.1",
    author="Paul Wambo",
    author_email="adjon081@uottawa.ca",
    description="Genomic Pathway Miner Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/poloarol/pathway-finder",
    packages=setuptools.find_packages(),
    classfiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
)