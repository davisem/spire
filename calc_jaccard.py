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
from collections import defaultdict
from abc import ABCMeta, abstractmethod

from hasher import HashFactory, Hasher
from kmerize import KMERPairFactory, KMERSingleFactory


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
		"""
		Where the magic happens
		:param list data: A list of integer hashes to recieve a second random hashing"""
		
		signature = []
		
		for hasher in self._hfuncs:
			#This is why it's called MinHash :D
			hashes = [hasher(x) for x in data]
			minhash = min(filter(lambda x: x != 0, hashes))
			signature.append(minhash)
		
		return np.array(signature)


class FullSignatureWorkflow(object):
	IKMER_FACTORY = KMERPairFactory
	SEED = 10
	
	def __init__(self, kmer_size, n_hash_functions):
		self._n_hash_functions = n_hash_functions
		self._kmer_size = kmer_size

	def flow(self, seq):
		#Set everything up
		hash_factory = HashFactory(Hasher, self._n_hash_functions, self.SEED)
		hashes = hash_factory.getHashes()
		m_hasher = MinHash(h_funcs=hashes)

		#Now calculate
		hashed_values = self.IKMER_FACTORY().create(seq, self._kmer_size)
		return m_hasher.calcMinHash(hashed_values)


	def getJaccard(self, sig1, sig2):
		"""Actually calculates the Jaccard index"""

		return float(sum(sig1 == sig2)) / float(self._n_hashes)

class FastSignatureWorkflow(object):
	
	IKMER_FACTORY = KMERSingleFactory

	def __init__(self, ikmer_factory=IKMER_FACTORY):
		self._ikmer_factory = ikmer_factory

	def flow(self):
		pass


def workflow():
	"""Just do something for sanity check"""
	
	#Some test sequences
	seq1 = 'AAAAATTTTTTTTTCCCCCCCCC'* 100
	seq2 = 'AAAAATTTTTTTTTCCCCCCCCC'* 99 + 'AAAAATTCCCTTTTCCCCCCCCC'
	
	#Make our random hashing functions
	seed = 10
	hash_factory = HashFactory(Hasher, 200, seed)
	
	hashes = hash_factory.getHashes()
	m_hasher = MinHash(h_funcs=hashes)

	#Make some kmers from the test sequences
	seq1_kmers = KMERPairFactory.create(seq1, 10)
	seq2_kmers = KMERPairFactory.create(seq2, 10)
	
	# Get signatures from the kmers
	signature_1 = m_hasher.calcMinHash(seq1_kmers)
	signature_2 = m_hasher.calcMinHash(seq2_kmers)

	#3lau

	print m_hasher.getJaccard(signature_1, signature_2)


def main():
	return workflow()

if __name__ == "__main__":
	main()





