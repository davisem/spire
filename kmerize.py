#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'kmerize'


import sys
import numpy as np
from collections import defaultdict
from abc import ABCMeta, abstractmethod


class IKMERFactory(object):

	def __init__(self):
		self.words = self.makeWords()

	def create(self, w, k):
		return self.getHashValues()


class KMERPairFactory(IKMERFactory):
	"""Class for enumerating kmers from a window, c/o Patrick"""

	@classmethod
	def makeWords(cls, w, k):
		pairs = {}
		for i in range (len(w) - k):
			for j in range (i + 1, len(w) - k):
				a = w[i:i + k]
				b = w[j:j + k]
				d = j - i
				ab = a + b
				if ab in pairs:
				    pairs[ab] = min(d, pairs[ab])
				else :
				    pairs[ab] = d

		return pairs

	def getHashValues(self):
		"""Formats the kmer_pair, distance dict into a list of inter values"""

		return [hash(keys+str(distances)) for keys, distances in self.words.iteritems()]


class KMERSingleFactory(IKMERFactory):

	@classmethod
	def makeWords(cls, w, k):
		
		kmers = defaultdict(list)
		
		for i in xrange(len(w) - k + 1):
			kmers[w[i: i + k]].append(i)

		return kmers

	def getHashValues(self, words):
		return [hash(keys) for keys in words.iterkeys()]
