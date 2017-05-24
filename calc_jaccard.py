#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = "wat wat"
__modname__ = 'calc_jaccard'


import sys
import numpy as np

from abc import ABCMeta, abstractmethod


class KMERPairFactory(object):
	"""Class for enumerating kmers from a window, c/o Patrick"""

	@staticmethod
	def create(w, k):
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


class IHasher(object):
	"""
	Base class for a hash function
	"""
	
	__metaclass__ = ABCMeta

	def __call__(self, x):
		"""
		Make the class callable
		:param x: the value to hash.
		:return int: The hash value of the input.
		"""
		
		return self.hashIt(x)

	@abstractmethod		
	def hashIt(self, x):
		"""
		To be implemented by child classes.
		"""
		pass


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
		return ((self._a * x) + self._b) % self._c


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
		
		randoms = np.random.RandomState(seed=self._seed).randint(self.MININT, self.MAXINT, self.n_hashes)
		
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


class MinHash(object):
	"""A class for implementing the minhash functionality"""


	def __init__(self, h_funcs):
		"""
		Init method for class
		:param list h_funcs: A list of classes implementing the Hasher interface"
		"""
		self._hfuncs = h_funcs
		self._n_hashes = len(self._hfuncs)


	def calcMinHash(self, data):
		"""Where the magic happens"""
		
		signature = []
		
		data = self.getData(data)
		
		for hasher in self._hfuncs:
			#This is why it's called MinHash :D
			hashes = [hasher(x) for x in data]
			minhash = min(filter(lambda x: x != 0, hashes))
			signature.append(minhash)
		
		return np.array(signature)

	@staticmethod
	def getData(kmer_pairs):
		"""Formats the kmer_pair, distance dict into a list of inter values"""

		return [hash(keys+str(distances)) for keys, distances in kmer_pairs.iteritems()]

	def getJaccard(self, sig1, sig2):
		"""Actually calculates the Jaccard index"""

		return float(sum(sig1 == sig2)) / float(self._n_hashes)


def workflow():
	"""Just do something for sanity check"""
	
	#Some test sequences
	seq1 = 'AAAAATTTTTTTTTCCCCCCCCC'
	seq2 = 'AAAAATTTTTTTTTCCCCCCCCC'
	
	#Make our random hashing functions
	seed = 10
	hash_factory = HashFactory(Hasher, 200, seed)
	
	hashes = hash_factory.getHashes()
	m_hasher = MinHash(h_funcs=hashes)

	#Make some kmers from the test sequences
	seq1_kmers = KMERPairFactory.create(seq1, 4)
	seq2_kmers = KMERPairFactory.create(seq2, 4)
	
	# Get signatures from the kmers
	signature_1 = m_hasher.calcMinHash(seq1_kmers)
	print signature_1
	signature_2 = m_hasher.calcMinHash(seq2_kmers)

	#3lau

	print m_hasher.getJaccard(signature_1, signature_2)


def main():
	return workflow()

if __name__ == "__main__":
	main()





