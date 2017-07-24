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
import math
import numpy as np
from collections import defaultdict
from abc import ABCMeta, abstractmethod, abstractproperty


class Read(object):
	"""Class representing a read"""

	__slots__ = ['seq', '_sketch']

	def __init__(self, seq):
		self.seq = seq
		self._sketch = None

	def setSketch(self, sketch):
		self._sketch = sketch

	@property
	def sketch(self):
		return self._sketch


class Kmer(object):
	"""Class representing a kmer"""
	
	__metaclass__ = ABCMeta

	@abstractproperty
	def integer_val(self):
		pass

	@abstractproperty
	def idx(self):
		pass

	@abstractproperty
	def __str__(self):
		pass


class KmerSingle(Kmer):
	"""Class representing a single kmer"""
	
	__slots__ = ['seq', 'idx', 'integer_val']
	
	def __init__(self, idx, seq):
		self.seq = seq
		self.integer_val = hash(seq)
		self.idx = idx

	def __str__(self):
		return self.seq

	def __eq__(self, other):
		return self.seq == other.seq

	def __repr__(self):
		return "<Kmer:{}, idx:{}>".format(self.seq, self.idx)

	def __index__(self):
		return self.idx

	def __gt__(self, other):
		return self.idx > other


class KmerPair(Kmer):
	"""Class representing a kmer pair"""

	__slots__ = ['kmer1', 'kmer2', 'distance', 'integer_val']

	def __init__(self, kmer1, kmer2):
		self.kmer1 = kmer1.seq
		self.kmer2 = kmer2.seq
		self.distance = kmer2.idx - kmer1.idx
		self.integer_val = hash(self.kmer1 + self.kmer2 + str(self.distance))

	@property
	def idx(self):
		pass

	def __str__(self):
		pass

	def __eq__(self, other):
		self.integer_val == other.integer_val

	def __gt__(self, other):
		return self.distance > other.distance


class Window(object):
	"""A class that holds logic for kmerizing a string"""
	
	KMER = KmerSingle

	def __init__(self, seq, kmer_size, n_buckets):
		"""
		Init method for class
		:param str seq: Input sequence
		:param int kmer_size: The size of kmers to enumerate
		:param n_buckets: Threshold for determining the distance metric if kmer-pairs are enumerated
		"""
		self.seq = seq
		self.buckets = self.getBuckets(n_buckets, len(seq))
		self.k = kmer_size


	def __getitem__(self, key):
		return self.kmers[key]

	def makeWords(self, word):
		"""
		Enumerates the window into kmers
		:param str word: The word to enumerate
		:param int k: The size of kmers to enumerate
		return list(KmerSingle): A list of enumerated KmerSingle instances
		"""
		return [self.KMER(self.buckets[i], word[i: i + self.k]) for i in xrange(len(word) - self.k + 1)]

	def __repr__(self):
		return "<{}:n_kmers:{}>".format(self.__class__.__name__, len(self.kmers))

	def getDownstreamKmers(self, kmer, kmers, distance):
		"""
		Gets all kmers downstream of a given kmer. Will use the first kmer found in the list
		:param Kmer kmer: The kmer object to enumerate downstream from
		"""
		for kmer in self.kmers[index:distance]:
			yield kmer

	def makeKmerPairs(self, word, distance):
		"""
		Enumerates kmer pairs
		:param str word: The string to enumerate
		:param int distance: A distance threshold to server as a bound on the kmer pair enumeration.
		:return list(KmerPair): A list of KmerPairs
		"""
	    

		kmerpairs = {}

		distance = min(len(self.seq), distance)
		
		for i in xrange(len(word) - self.k + 1):
			
			for j in xrange(i+1, distance):
				
				first_kmer = KmerSingle(self.buckets[i], word[i:i + self.k])
				second_kmer = KmerSingle(self.buckets[j], word[j:j + self.k])
				pair = KmerPair(first_kmer, second_kmer)
				
				if pair.integer_val in kmerpairs:
					eval_pair = kmerpairs[pair.integer_val]
					kmerpairs[pair.integer_val] = min(eval_pair, pair)
				
				else :
					kmerpairs[pair.integer_val] = pair
		
		return kmerpairs.values()

	@staticmethod
	def getBuckets(n_buckets, len_seq):
		"""
		Sets up a bucketing scheme for which kmer distances can be measured:
		:param int n_buckets: The number of buckets to consider for the sequence
		:param int len_seq: How long the sequence is to divide into buckets
		:return dict: A mapping of index to bucket value
		"""
		items_per_bucket = math.ceil(float(len_seq) / float(n_buckets))
		running_bucket = range(n_buckets)
		current_bucket = running_bucket[0]
		buckets = {}
		bucket_count = 0
		
		for i in range(len_seq):
			
			if bucket_count < items_per_bucket:
				
				buckets[i] = running_bucket[current_bucket]
				bucket_count += 1
			
			else:
				current_bucket += 1
				bucket_count = 0
				buckets[i] = running_bucket[current_bucket]
		
		return buckets
