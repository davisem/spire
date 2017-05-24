#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'fastq_filter_tests.py'


import unittest

from calc_jaccard import MinHash, KMERize


class TestMinHash(unittest.TestCase):


	def test_minhash(self):
		pass



class TestKMERize(unittest.TestCase):

	def test_kmerize_keys(self):
		test = 'AAAATTTT'
		expect_keys = ['AAT', 'ATT', 'AAA', 'TTT']
		kmers = KMERize.create(test, 3)
		self.assertListEqual(kmers.keys(), expect_keys)

	def test_kmerize_values(self):
		test = 'AAAATTTT'
		expect_values = [0, 1, 2, 3, 4, 5]
		kmers = KMERize.create(test, 3)
		obs_values =[]
		for positions in kmers.values():
			obs_values.extend(positions)
		self.assertListEqual(sorted(obs_values), expect_values)

