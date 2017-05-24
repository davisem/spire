#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'hasher'

import numpy as np

class Hasher(IHasher):
	"""Some hashing function"""
	
	__slots__ = ['a', 'b']

	def __init__(self, a, b):
		self._a = a
		self._b = b
		self._c = max(a, b) - 100

	def hashIt(self, x):
		"""
		Hashes the input.
		:param int x: The input value to hash.
		:return int: The hash value.
		"""
		return (x + self._b) % self._c


class HashFactory(object):
	"""
	h(x) = (ax + b) % c
	such that a and b are random ints < x. 
	"""

	MININT = np.iinfo(np.int32).min
	MAXINT = np.iinfo(np.int32).max

	def __init__(self, ihasher, n_hashes, seed):
		"""
		Init method for class
		:param IHasher ihasher: A class implementing the ihasher interface
		:param int n_hashes: How many random hashes to make
		:param int seed: Seed for random number generation
		"""
		
		self.ihasher = ihasher
		self.n_hashes = n_hashes
		self._seed =  seed

	def getRandomInts(self):
		"""
		Makes a pair of random integers. 
		Seed ensures the randomness is deterministic.
		"""
		
		randoms = np.random.RandomState(seed=self._seed).randint(self.MININT, self.MAXINT, self.n_hashes * 2)
		
		for i in range(len(randoms) - 2):
			a, b = randoms[i: i + 2]
			yield a, b

	def getHashes(self):
		"""
		Returns instances of random hashing classes.
		:return list(IHasher): a list of IHasher instances.
		"""
		h_funcs = []
		for a, b in self.getRandomInts():
			 h_funcs.append(self.ihasher(a, b))
		
		return h_funcs