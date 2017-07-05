#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = "min_hash.py"

import sys
import numpy as np
import operator
from collections import defaultdict, Counter
from abc import ABCMeta, abstractmethod

from hasher import HashFactory
from kmer import Kmer, Window, Read

import pyximport; pyximport.install()
from cython.chash import Hasher


class MinHash(object):
	"""A class for implementing the minhash functionality"""

	def __init__(self, n_hash_functions):
		"""
		Init method for class
		:param list h_funcs: A list of classes implementing the Hasher interface"
		"""
		self._n_hash_functions = n_hash_functions
		self._hash_permutations = HashFactory(Hasher, n_hash_functions).getHashes()

	@staticmethod
	def getMinimumHash(hashes_kmer):
		return min(hashes_kmer, key=operator.itemgetter(0))

	def calcSketch(self, kmers):
		"""
		Calculates a minhash Sketch given a set of kmers.
		:param list(Kmer) kmers: a list of Kmer instances
		:return list(int) Sketch: A list representing a minhash Sketch
		"""
		
		sketch = []

		for hasher in self._hash_permutations:
			
			hashes = [(hasher.hashIt(kmer.integer_val), kmer) for kmer in kmers]
			min_hash = self.getMinimumHash(hashes)
			
			sketch.append(min_hash)
		
		return sketch

	def approximateJaccard(self, sig1, sig2):
		"""
		Calculate an approximate Jaccard index based on minhash Sketchs
		TODO: Assert that sig1 and sig2 were calculated with the same hash functions, otherwise
		the approimation is invalid.
		:param list(int) sig1:
		:param sig2
		"""
		
		sig1 = np.array(sig1)
		sig2 = np.array(sig2)
		return float(sum(sig1 == sig2)) / float(self._n_hash_functions)


class SketchDirector(object):
	"""Flow control for generating minhash Sketchs for a given sequence"""

	__metaclass__ = ABCMeta

	def __init__(self, word_size, n_hashing_functions, hasher):
		"""
		Init method for class
		:param int word_size: The kmer size to use
		:param n_hashing_functions: The number of hashing functions to use.
		"""
		self.word_size = word_size
		self._n_hashing_functions = n_hashing_functions
		self._min_hash = hasher(n_hashing_functions)
		self._Sketch_cache = {x:defaultdict(list) for x in xrange(n_hashing_functions)}

	def run(self, seq):
		Sketch = self.getSketch(seq)
		self.storeSketch(seq, Sketch)

	def storeSketch(self, seq, Sketch):
		for i, min_mer in enumerate(Sketch):
			self._Sketch_cache[i][min_mer].append(seq)

	def querryRead(self, Sketch):
		for i, min_mer in enumerate(Sketch):
			yield self._Sketch_cache[i][min_mer]

	def getSimilarReads(self, read):
		
		similar_reads = []
		for q_read_list in self.querryRead(read.sketch):
			for q_read in q_read_list:
				if q_read != read.seq:
					similar_reads.append(q_read)
		
		return Counter(similar_reads)


class GreedyPairSketchDirector(SketchDirector):


	def __init__(self, word_size, n_hashing_functions, hasher):
		super(GreedyPairSketchDirector, self).__init__(word_size, n_hashing_functions, hasher)
	
	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		
		window = Window(seq, self.word_size, 2)
		
		primary_Sketch = self._min_hash.calcSketch(window.kmers)
		kmer_pairs = []
		
		for min_hash in primary_Sketch:
			kmer_pairs.extend(window.makeKmerPairs(min_hash[1]))
		
		return [x[0] for x in self._min_hash.calcSketch(kmer_pairs)]


class SingleSketchDirector(SketchDirector):
	
	def __init__(self, word_size, n_hashing_functions, hasher):
		super(SingleSketchDirector, self).__init__(word_size, n_hashing_functions, hasher)

	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		window = Window(seq, self.word_size, 1)
		return [x[0] for x in self._min_hash.calcSketch(window.kmers)]


class ExhaustivePairSketchDirector(SketchDirector):

	def __init__(self, word_size, n_hashing_functions, hasher):
		super(ExhaustivePairSketchDirector, self).__init__(word_size, n_hashing_functions, hasher)


	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		
		window = Window(seq, self.word_size, 2)
		
	
		pairs = [pair for kmer in window.kmers for pair in window.makeKmerPairs(kmer)]
		return [x[0] for x in self._min_hash.calcSketch(pairs)]
