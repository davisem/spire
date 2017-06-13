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
import operator
from collections import defaultdict
from abc import ABCMeta, abstractmethod

from hasher import HashFactory, Hasher
from kmerizer import KMERPairFactory, KMERSingleFactory
from min_hash import MinHash


class SignatureCalculator(object):
	"""Class for calcluating a unique minhash signature from a given word, or window of a sequence"""

	__metaclass__ = ABCMeta

	def __init__(self, kmer_size, hasher, n_hash_functions, kmer_factory):
		self._n_hash_functions = n_hash_functions
		self._hash_permutations = HashFactory(hasher, self._n_hash_functions).getHashes()
		self._min_hasher = MinHash()
		
		self._kmer_factory = kmer_factory(kmer_size)
	
	def calcJaccard(self, sig1, sig2):
		"""Actually calculates the Jaccard index"""
		sig1 = np.array(sig1)
		sig2 = np.array(sig2)
		return float(sum(sig1 == sig2)) / float(self._n_hash_functions)

	@abstractmethod
	def calculate(self, window):
		pass


class PairwiseDistanceSignatureCalculator(SignatureCalculator):
	"""Calclulates a minhash signature using kmer pairs, and a distance value. SLOW"""
	
	IKMER_FACTORY = KMERPairFactory

	def __init__(self, kmer_size, hasher, n_hash_functions, kmer_factory=IKMER_FACTORY):
		super(PairwiseDistanceSignatureCalculator, self).__init__(kmer_size, hasher, n_hash_functions, kmer_factory)
	
	def calculate(self, window):
		
		signature = []

		kmers = self._kmer_factory.create(window)

		for hasher in self._hash_permutations:
			min_hash = self._min_hasher.getMinimumHash(kmers, hasher)
			signature.append(min_hash)
		
		return signature


class FastSignatureCalculator(SignatureCalculator):
	"""Calculates a signature by looking up the kmer-pair and distance value only from the minhash value"""

	IKMER_FACTORY = KMERSingleFactory

	def __init__(self, kmer_size, hasher, n_hash_functions):
		super(FastSignatureCalculator, self).__init__(kmer_size, hasher, n_hash_functions, FastSignatureCalculator.IKMER_FACTORY)

	def calculate(self, window):
		
		signature = []
		
		kmers = self._kmer_factory(window, self._kmer_size)
		
		for hasher in self._hash_permutations:
			
			min_hash = self._min_hasher.getMinimumHash(hasher, kmers)
			signature.append(min_hash)
		
		return signature

	def factorDistances(self):
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





