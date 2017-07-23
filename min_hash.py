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
from hashers.chash import Hasher


class MinHash(object):
	"""A class for implementing the minhash functionality"""

	def __init__(self, n_hash_functions):
		"""
		Init method for class
		:param list h_funcs: A list of classes implementing the Hasher interface"
		"""
		self._n_hash_functions = n_hash_functions
		self._hash_permutations = HashFactory(Hasher, n_hash_functions).getHashes()

	def calcSketch(self, kmers):
		"""
		Calculates a minhash Sketch given a set of kmers.
		:param list(Kmer) kmers: a list of Kmer instances
		:return list(int) Sketch: A list representing a minhash Sketch
		"""
		
		sketch = []
		i = 0
		print("Calculating {} hashes across {} kmers".format(len(self._hash_permutations), len(kmers)))
		
		for hasher in self._hash_permutations:
			if i % 100 == 0:
				print ("Calculated for {} hash functions so far".format(i))
			
			sketch.append(min([hasher.hashIt(kmer.integer_val) for kmer in kmers]))
			
			i += 1
		
		return sketch
