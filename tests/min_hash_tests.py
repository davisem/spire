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
import cProfile
from pstats import Stats

from long_read_aligner.min_hash import MinHash, LoggingMinHash, SignatureDirector
from long_read_aligner.kmer import Kmer, Window
from long_read_aligner.utils.string_utils import CreateWindow, CreateKMERPairs


class MinHashTests(unittest.TestCase):

	def test_calcSignature(self):
		"""Test that the signature lenght is equal to the number of hashing fucntions"""
		n_hash_permutations = 10
		l_min_hash = LoggingMinHash(n_hash_permutations)

		kmers = ['AA', 'TT', 'CC', 'GG']
		kmers = [Kmer(i, x) for i, x in enumerate(kmers)]

		signature = l_min_hash.calcSignature(kmers)
		self.assertEqual(len(signature), n_hash_permutations)

	def test_getMinimumHash(self):
		"""Test that we can get the minimum hash"""
		hashes = [0, 1, 2, 3, 4, 5]
		min_hash = MinHash.getMinimumHash(hashes)
		self.assertEquals(1, min_hash)

	def test_approximateJaccardGood(self):
		"""Test that the approximate jaccard is 1 at max identity"""
		min_hash = MinHash(3)
		sig1 = [20, 20, 10]
		sig2 = [20, 20, 10]
		jaccard = min_hash.approximateJaccard(sig1, sig2)
		self.assertEquals(1, jaccard)

	def test_approximateJaccardBad(self):
		min_hash = MinHash(3)
		sig1 = [100, 200, 30]
		sig2 = [20, 20, 10]
		jaccard = min_hash.approximateJaccard(sig1, sig2)
		self.assertEquals(0, jaccard)


class SignatureDirectorTests(unittest.TestCase):

	def setUp(self):
		"""Enable profiling"""
		self.pr = cProfile.Profile()
		self.pr.enable()
		print "\n<<<---"

	def tearDown(self):
		"""Report profiling results"""
		p = Stats(self.pr)
		p.strip_dirs()
		p.sort_stats('cumtime')
		p.print_stats()
		print "\n--->>>"

	def test_getSignature(self):
		""""""
		test_seq = "AAATTACCCGGG"
		test_seq2 = "AAATTTCCCGGG"
		director = SignatureDirector(3, 100)
		
		sig1 = director.getSignature(test_seq)
		sig2 = director.getSignature(test_seq2)

		jaccard = director._min_hash.approximateJaccard(sig1, sig2)
		print jaccard
		self.assertTrue(0 < jaccard < 1)

	def test_time(self):
		"""Kind of a stress test to get some time benchmarks. Here we generate a signature for a 1kb window"""
		test_seq = CreateWindow(1000)
		director = SignatureDirector(8, 100)
		director.getSignature(test_seq)



