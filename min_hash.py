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
from collections import defaultdict
from abc import ABCMeta, abstractmethod

from hasher import HashFactory, Hasher
from kmer import Kmer, Window


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
	def getMinimumHash(hashes):
		return min(filter(lambda x: x != 0, hashes))

	def calcSignature(self, kmers):
		"""
		Calculates a minhash signature given a set of kmers.
		:param list(Kmer) kmers: a list of Kmer instances
		:return list(int) signature: A list representing a minhash signature
		"""
		
		signature = []

		for hasher in self._hash_permutations:
			hashes = [hasher(x) for x in kmers]
			min_hash = self.getMinimumHash(hashes)
			signature.append(min_hash)
		
		return signature

	def approximateJaccard(self, sig1, sig2):
		"""
		Calculate an approximate Jaccard index based on minhash signatures
		TODO: Assert that sig1 and sig2 were calculated with the same hash functions, otherwise
		the approimation is invalid.
		:param list(int) sig1:
		:param sig2
		"""
		
		sig1 = np.array(sig1)
		sig2 = np.array(sig2)
		return float(sum(sig1 == sig2)) / float(self._n_hash_functions)


class LoggingMinHash(MinHash):
	"""
	Class for calculating a min hash signature. Implements logging capabilities to back-lookup original input values
	prior to hashing. This removes the need to calculate all kmer pairs and distances exhaustively, which is a 
	quadratic operation. Alot of book keeping going on here.
	"""

	def __init__(self, n_hash_functions):
		super(LoggingMinHash, self).__init__(n_hash_functions)
		self.hash_to_kmer = {}

	def setHashToKmer(self, min_hash_to_int):
		"""Updates the hash to int map of the instance"""
		self.hash_to_kmer.update(min_hash_to_int)

	def getKmer(self, min_hash):
		"""
		Gets the kmer string given the minhash value.
		"""
		return self.hash_to_kmer[min_hash]

	@staticmethod
	def getMinimumHash(hashes_kmer):
		return min(hashes_kmer, key=operator.itemgetter(0))

	def calcSignature(self, kmers):
		"""
		Calculates a minhash signature given a set of kmers.
		:param list(Kmer) kmers: A list of Kmers instances
		"""

		primary_signature = []

		for hasher in self._hash_permutations:

			hashes = [(hasher.hashIt(x.integer_val), x) for x in kmers if x.integer_val!=0]
			
			min_hash = self.getMinimumHash(hashes)
			
			primary_signature.append(min_hash[0])
			
			self.hash_to_kmer[min_hash[0]] = min_hash[1]

		return primary_signature


class SignatureDirector(object):
	"""Flow control for generating minhash signatures for a given sequence"""

	def __init__(self, word_size, n_hashing_functions):
		"""
		Init method for class
		:param int word_size: The kmer size to use
		:param int n_hashing_functions: The number of hashing functions to use.
		"""
		self.word_size = word_size
		self._min_hash = LoggingMinHash(n_hashing_functions)

	def getSignature(self, seq):
		"""
		Gets the minhash signature for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash signature
		"""
		
		window = Window(seq, self.word_size)
		primary_signature = self._min_hash.calcSignature(window.kmers)
		kmer_pairs = []
		
		for min_hash in primary_signature:
			min_kmer = self._min_hash.getKmer(min_hash)
			kmer_pairs.extend(window.makeKmerPairs(min_kmer))

		return self._min_hash.calcSignature(kmer_pairs)




		




