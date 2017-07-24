#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'kmer_tests.py'


import unittest

from kmer import KmerSingle, KmerPair, Window


class TestKmerSingle(unittest.TestCase):

	def setUp(self):
		self.kmer = KmerSingle(0, "AAAT")

	def test_Str(self):
		"""Test that we can get the kmer in string format"""
		self.assertEqual(self.kmer.seq, str(self.kmer))

	def test_GreaterThan(self):
		"""Test that the greater than functionality works according to the index"""
		greater_kmer = KmerSingle(2, "AAA")
		self.assertFalse(greater_kmer > self.kmer)


class TestWindow(unittest.TestCase):

	def setUp(self):
		self.window = Window("AAAT", 2, 2)

	def test_makeWords(self):
		"""Test that kmers are generated correctly"""
		kmer_values = ['AA', 'AA', 'AT']
		expect = [KmerSingle(i, x) for i, x in enumerate(kmer_values)]
		words = self.window.makeWords("AAAT")

	def test_makePairs(self):
		"""Test that the pairs are correctly generated"""
		k1 = KmerSingle(0, 'AA')
		k2 = KmerSingle(1, 'AC')
		k3 = KmerSingle(2, 'CT')
		
		p1 = KmerPair(k1, k2)
		p2 = KmerPair(k1, k3)
		p3 = KmerPair(k2, k3)

		expect = [p1, p2, p3]
		window = Window("AACT", 2, 2)
		words = window.makeKmerPairs("AACT", 3)
		self.assertTrue(len(expect), len(words))

