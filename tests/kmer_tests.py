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

from long_read_aligner.kmer import Kmer, Window


class TestKmer(unittest.TestCase):

	def setUp(self):
		self.kmer = Kmer(0, "AAAT")

	def test_Str(self):
		"""Test that we can get the kmer in string format"""
		self.assertEqual(self.kmer.seq, str(self.kmer))

	def test_GreaterThan(self):
		"""Test that the greater than functionality works according to the index"""
		greater_kmer = Kmer(2, "AAA")
		self.assertFalse(greater_kmer > self.kmer)


class TestWindow(unittest.TestCase):

	def setUp(self):
		self.window = Window("AAATTTCCCGGG", 3)

	def test_makeWords(self):
		"""Test that kmers are generated correctly"""
		kmer_values = ['AA', 'AA', 'AT']
		expect = [Kmer(i, x) for i, x in enumerate(kmer_values)]
		words = Window.makeWords('AAAT', 2)
		self.assertListEqual(words, expect)

	def test_getDownstreamKmers(self):
		"""Test that we can get the next downstream kmer from the string"""
		kmer = Kmer(0, "AAA")
		for i, kmer in enumerate(self.window.getDownstreamKmers(kmer)):
			self.assertEqual(kmer, self.window[i+1])

	def test_makeKmerPairs(self):
		"""Test that we can build kmer pairs correctly."""
		kmer = Kmer(0, "AAT")
		kmer_pairs = self.window.makeKmerPairs(kmer)
		self.assertEquals(kmer_pairs[0].seq, "AATATT2")
