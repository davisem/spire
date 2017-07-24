#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'min_hash_tests.py'

import unittest

from pstats import Stats

from kmer import Read
from kmer import KmerSingle, Window
from min_hash import MinHash
from utils.string_utils import CreateWindow, CreateKMERPairs, PerturbWindow


class MinHashTests(unittest.TestCase):

	def test_calcSketchSingle(self):
		"""Test that the Sketch lenght is equal to the number of hashing fucntions"""
		n_hash_permutations = 10
		l_min_hash = MinHash(n_hash_permutations)

		kmers = ['AA', 'TT', 'CC', 'GG']
		kmers = [KmerSingle(i, x) for i, x in enumerate(kmers)]

		Sketch = l_min_hash.calcSketch(kmers)
		self.assertEqual(len(Sketch), n_hash_permutations)
