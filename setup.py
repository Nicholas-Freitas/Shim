from setuptools import setup, find_packages

setup(
    name='shim',
    version='0.0.1',    
    description='A utility for fixing and converting molecular structures',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Nicholas-Freitas/shim',
    author='Nicholas Freitas',
    author_email='nicholas.freitas@ucsf.edu',
    license='MIT',
    packages=find_packages(),
    entry_points={},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=[
    'rdkit',
    'biopython',
    'termol>=0.1.7' #this is causing issues due to rdkit/rdkit-pypi reqs.
    ],
)
