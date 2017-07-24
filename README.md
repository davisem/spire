[![CircleCI](https://circleci.com/gh/davisem/spire.svg?style=shield)](https://circleci.com/gh/davisem/spire/master)

Spire
=========================
An experimental, long read genome-assembler that leverages k-mer pairs along with distance information to map error-prone reads with high fidelity. This application is very much a work in progress, and consists of a simple command line interface for robustly approximating a Jaccard index between two strings. More to come...

Prerequisites
-------------
  - Python 2.7
  - Cython
  - Numpy

Quickstart
==========

Clone Github Repository
-----------------------
```
$ git clone https://github.com/davisem/spire.git
$ cd spire
```

Setup
------------------
```
$ virtualenv venv
$ source venv/bin/activate
$ pip install cython numpy
$ python setup.py build_ext --inplace
```

Usage
-----
```
$ python spire.py -q ATCGGATCGAATCG -r ATCGATCGATCGATCG -k 3 -m deep -n 300
```

Options
----------------
Option | Description
--------- | -----------
help | `Display help message`
reference | `A reference string to search against`
query_string| `A string to compare to the reference`
n_hashing_functions | `Number of random hashing functions to use`
word_size | `Word size with which input strings will be enumerated`
mode | `[fast, deep] Whether to enumerate kmers (fast) or kmer-pairs (deep)`

Contributors
============

Original concept by Hamidreza Chitsaz. Implemented by Chris Dean, Eric Davis, and Patrick Hagarty

Contact
=======
Feel free to send me an email at emdavis48 AT gmail DOT com