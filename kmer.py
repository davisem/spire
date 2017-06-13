#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = 'kmer.py'


import sys

class Kmer(object):
	
	__slots__ = ['seq', '_idx', 'integer_val']
	
	def __init__(self, idx, seq):
		self.seq = seq
		self.integer_val = hash(seq)
		self._idx = idx

	def __str__(self):
		return self.seq

	def __eq__(self, other):
		return self.seq == other.seq

	def __repr__(self):
		return "<Kmer:{}, idx:{}>".format(self.seq, self._idx)

	def __index__(self):
		return self._idx

	def __gt__(self, other):
		return self._idx > other


class Window(object):
	
	KMER = Kmer

	def __init__(self, seq, kmer_size):
		self.seq = seq
		self.kmers = self.makeWords(seq, kmer_size)
		self.kmer_pairs = None

	def __getitem__(self, key):
		return self.kmers[key]

	@classmethod
	def makeWords(cls, w, k):
		return [cls.KMER(i, w[i: i + k]) for i in xrange(len(w) - k + 1)]

	def __repr__(self):
		return "<{}:n_kmers:{}>".format(self.__class__.__name__, len(self.kmers))

	def getDownstreamKmers(self, kmer):
		index = self.kmers.index(kmer) + 1
		for kmer in self.kmers[index:]:
			yield kmer

	def makeKmerPairs(self, kmer):
		kmer_idx = kmer._idx
		return [Kmer(kmer_idx, ''.join([kmer.seq, next_kmer.seq, str(next_kmer._idx - kmer_idx)])) for next_kmer in self.getDownstreamKmers(kmer)]


