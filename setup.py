from setuptools import setup, find_packages

setup(
    name="metaCH",
    version="0.1",
    packages=find_packages(where="metaCH"),
    package_dir={"": "metaCH"},
    install_requires=[],
    author="Marzieh Haghighi",
    description="MetaCH: An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples",
)
