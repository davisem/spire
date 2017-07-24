#!/usr/bin/env python

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__modname__ = 'spire.py'


import argparse
import sys
import numpy as np
from min_hash import MinHash
from sketch_calculators import SingleSketchCalculator, ExhaustivePairSketchCalculator
from utils.string_utils import PerturbWindow


def parse_cmdline_params(cmdline_params):
    
    info = "Calculate a Jaccard index for a query sequence against a reference"
    
    parser = argparse.ArgumentParser(description=info)
    
    parser.add_argument('-q', '--query_string', type=str, required=True,
                        help='Provide a query string')

    parser.add_argument('-r', '--reference_sequence', type=str, required=True,
    					help='Please provide a reference string')

    parser.add_argument('-n', '--n_hashing_functions', type=int, default=100, 
    					help='Number of hashing f(x) to considering when calculating a signature')

    parser.add_argument('-k', '--word_size', type=int, help="kmer_size")

    parser.add_argument('-m', '--mode', type=str, choices=['greedy', 'fast', 'deep'], 
    					default='fast', help="The algorithm to generate the kmer strings for hashing")

    return parser.parse_args(cmdline_params)


class Spire(object):
	
	"""Spire application"""

	CALCULATORS = [SingleSketchCalculator, ExhaustivePairSketchCalculator]

	def __init__(self, word_size, n_hash_functions, mode):
		"""
		Init method for class
		:param int word_size: The size of words to kmerize
		:param int n_hash_fucntions: Number of random hash functions to use for the min_hash signature
		:param str mode: ["fast", "deep"] Whether to enumerate single kmers or kmer-pairs
		"""
		self._word_size = word_size
		self._n_hash_functions = n_hash_functions
		self._calculator = self._setMode(mode)
		self.ref_signature = None

	def _setMode(self, mode):
		"""
		Dispatches the correct mode
		:param str mode: ["fast", "deep"] Make single kmers or kmer pairs.
		"""
		for calculator in self.CALCULATORS:
			
			if calculator.canCalculate(mode):
			
				return calculator(self._word_size, self._n_hash_functions)

		raise ValueError("No calculator was found to hande mode: \"{}\"".format(mode))


	def loadRefSignature(self, ref):
		"""
		Calculate and store the reference signature
		:param str ref: The reference sequence
		"""
		self.ref_signature = self._calculator.getSketch(ref)

	def getJaccard(self, query_sequence):
		"""
		Get the jaccard similarity for a query sequence
		:param str query_sequence: The sequence to calculate the jaccard against the reference
		:return float: The jaccard
		"""
		query_signature = self._calculator.getSketch(query_sequence)
		return self.approximateJaccard(self.ref_signature, query_signature)

	def approximateJaccard(self, sig1, sig2):
		"""
		Calculate an approximate Jaccard index based on minhash Sketchs
		TODO: Assert that sig1 and sig2 were calculated with the same hash functions, otherwise
		the approimation is invalid.
		:param list(int) sig1: A list representation of a min_hash signature
		:param list(int) sig2: A list representation of a min_hash signature
		:return float: The jaccard
		"""
		sig1 = np.array(sig1)
		sig2 = np.array(sig2)
		return (float(sum(sig1 == sig2)) / float(self._n_hash_functions))


if __name__ == '__main__':

	modes = [SingleSketchCalculator, ExhaustivePairSketchCalculator]
	opts= parse_cmdline_params(sys.argv[1:])

	spire = Spire(opts.word_size, opts.n_hashing_functions, opts.mode)
	spire.loadRefSignature(opts.reference_sequence)

	jaccard = spire.getJaccard(opts.query_string)
	sys.stdout.write("Spire Approximated Jaccard: {}\n".format(jaccard))
