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

Pipeline Options
----------------
Option | Description
--------- | -----------
help | `Display help message.`
reference | `A reference string to search against`
querry_string| `A string to compare to the reference`
n_hashing_functions | `Number of random hashing functions to use`
word_size | `Word size with which input strings will be enumerated`
mode | `[Fast, Deep] Whether to enumerate kmers or kmer-pairs`


Tools and Core Libraries
========================
  - [Blasr](https://github.com/PacificBiosciences/blasr)
    - PacBio Long Read Alignerg
  - [Magic-Blast](https://www.ncbi.nlm.nih.gov/news/09-22-2016-magic-BLAST/)
    - Whole genome aligner
  - [Canu](https://github.com/marbl/canu)
    - Single molecule sequence assembler
  - [MHAP](https://github.com/marbl/MHAP)
    - Local sensitivity hashing

Contributors
============

Original concept by Hamidreza Chitsaz. Implemented by Chris Dean, Eric Davis, and Patrick Hagarty

Contact
=======
Feel free to send me an email at emdavis48 AT gmail DOT com