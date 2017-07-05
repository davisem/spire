#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'hasher.py'

import numpy as np



class HashFactory(object):
	"""
	Class for generating perumations of a hashing object.
	h(x) = (ax + b) % c
	such that a and b are random ints < x. 
	"""

	MININT = np.iinfo(np.int32).min
	MAXINT = np.iinfo(np.int32).max
	SEED = 10

	def __init__(self, ihasher, n_hashes):
		"""
		Init method for class
		:param IHasher ihasher: A class implementing the ihasher interface
		:param int n_hashes: How many random hashes to make
		:param int seed: Seed for random number generation
		"""
		
		self.ihasher = ihasher
		self.n_hashes = n_hashes

	def getRandomInts(self):
		"""
		Makes a pair of random integers. 
		Seed ensures the randomness is deterministic.
		"""
		
		randoms = np.random.RandomState(seed=self.SEED).randint(self.MININT, self.MAXINT, self.n_hashes + 1)
		
		for i in xrange(self.n_hashes):
			random_ints = randoms[i: i + 2]
			yield random_ints

	def getHashes(self):
		"""
		Returns instances of random hashing classes.
		:return list(IHasher): a list of IHasher instances.
		"""
		h_funcs = []
		for a, b in self.getRandomInts():
			 h_funcs.append(self.ihasher(a, b))
		
		return h_funcs
