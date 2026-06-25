#!/usr/bin/env python

import io
from setuptools import setup, find_packages

long_description=io.open("README.md",encoding='utf-8').read()

setup(name='stabgraph',
      version='0.1.5',
      description='Transforms stabilizer state into graph state',
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Science/Research",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3 :: Only",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Programming Language :: Python :: 3.11",
          "Programming Language :: Python :: 3.12",
          "Programming Language :: Python :: 3.13",
          "License :: OSI Approved :: Apache Software License",
          "Operating System :: OS Independent",
          "Topic :: Scientific/Engineering :: Physics",
      ],
      keywords='stabilizer graph local clifford',
      author='David Amaro',
      author_email='davamaro@hotmail.es',
      license='Apache 2',
      url='https://github.com/davamaro/stabgraph',
      project_urls={
          'Source': 'https://github.com/davamaro/stabgraph',
          'Issues': 'https://github.com/davamaro/stabgraph/issues',
      },
      packages=find_packages(),
      python_requires='>=3.9',
      install_requires=[
          'numpy',
          "galois>=0.4.11; python_version < '3.13'",
      ],
      extras_require={
          'accel': ["galois>=0.4.11; python_version < '3.13'"],
          'test': ['pytest'],
      },
      include_package_data=True)
