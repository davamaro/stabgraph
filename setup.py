#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
      long_description=fh.read()

setup(name='stabgraph',
      version='0.1',
      description='Transforms stabilizer state into graph state',
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "License :: OSI Approved :: Apache Software License",
          "Operating System :: OS Independent",
      ],
      keywords='stabilizer graph local clifford',
      author='David Amaro',
      author_email='davamaro@hotmail.es',
      license='Apache 2',
      packages=find_packages(),
      install_requires=[
          'numpy',
      ],
      include_package_data=True)
