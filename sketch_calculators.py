__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = "sketch_calculators.py"


import sys
import numpy as np
import operator
from collections import defaultdict, Counter
from abc import ABCMeta, abstractmethod

from hasher import HashFactory
from kmer import Kmer, Window, Read
from min_hash import MinHash

import pyximport; pyximport.install()
from hashers.chash import Hasher



class SketchCalculator(object):
	"""Flow control for generating minhash Sketchs for a given sequence"""

	MIN_HASH = MinHash
	N_BUCKETS = 2

	__metaclass__ = ABCMeta

	def __init__(self, word_size, n_hashing_functions):
		"""
		Init method for class
		:param int word_size: The kmer size to use
		:param n_hashing_functions: The number of hashing functions to use.
		"""
		self.word_size = self._setWordSize(word_size)
		self._n_hashing_functions = n_hashing_functions
		self._min_hash = self.MIN_HASH(n_hashing_functions)

	@abstractmethod
	def getSketch(self):
		"""Calculates a min_hash sketch"""
		pass

	@abstractmethod
	def _setWordSize(self, word_size):
		"""Set the word size from which kmers are enumerated"""
		pass


class SingleSketchCalculator(SketchCalculator):
	"""Calculates a min_hash sketch using single kmers"""

	MODE_KEY = "fast"

	def _setWordSize(self, word_size):
		"""
		Setter for the words size
		:param int word_size: The word size
		:return int: The word size to set
		"""
		return 2 * word_size

	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		window = Window(seq, self.word_size, self.N_BUCKETS)
		kmers = window.makeWords(seq)
		return self._min_hash.calcSketch(kmers)

	@classmethod
	def canCalculate(cls, mode):
		return mode == cls.MODE_KEY


class ExhaustivePairSketchCalculator(SketchCalculator):
	"""Calculate the sketch using kmer pairs plus distances"""
	
	MODE_KEY = "deep"
	PAIR_DISTANCE = 50

	def getSketch(self, seq):
		"""
		Gets the minhash Sketch for a given sequence
		:param str seq: A DNA sequence
		:return list(int): A list representing the minhash Sketch
		"""
		
		window = Window(seq, self.word_size, self.N_BUCKETS)
		pairs = window.makeKmerPairs(seq, self.PAIR_DISTANCE)
		return self._min_hash.calcSketch(pairs)

	def _setWordSize(self, word_size):
		"""
		Setter for the words size
		:param int word_size: The word size
		:return int: The word size to set
		"""
		return word_size

	@classmethod
	def canCalculate(cls, mode):
		return mode == cls.MODE_KEY
