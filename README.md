[![CircleCI](https://circleci.com/gh/davisem/spire.svg?style=shield)](https://circleci.com/gh/davisem/spire/master)

Spire
=========================
An experimental, long read genome-assembler that leverages k-mer pairs along with distance information to map highly error-prone reads from long read sequencing platforms with high fidelity.

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
$ cd long_read_aligner
```

Setup
------------------
```
$ python setup.py
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
mode | `[Fast, Deep] Whether to enumerate kmers or kmer-pairs`

Contributors
============

Original concept by Hamidreza Chitsaz. Implemented by Chris Dean, Eric Davis, and Patrick Hagarty

Contact
=======
Feel free to send me an email at emdavis48 AT gmail DOT com