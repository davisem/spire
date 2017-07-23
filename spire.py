

__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = "to the moon"

import argparse
import sys
import numpy as np
from min_hash import MinHash
from sketch_calculators import GreedyPairSketchCalculator, SingleSketchCalculator, ExhaustivePairSketchCalculator
from utils.string_utils import PerturbWindow

def parse_cmdline_params(cmdline_params):
    
    info = "Calculate a Jaccard index for a query sequence against a reference"
    
    parser = argparse.ArgumentParser(description=info)
    
    parser.add_argument('-q', '--query_string', type=str, required=True,
                        help='Provide a query string')

    parser.add_argument('-r', '--reference_sequence', type=str, required=True,
    					help='Please provide a reference string')

    parser.add_argument('-p', '--perturb_input', type=float, 
    					help='probabilty of a sequencing error at any given base')

    parser.add_argument('-n', '--n_hashing_functions', type=int, default=100, 
    					help='Number of hashing f(x) to considering when calculating a signature')

    parser.add_argument('-k', '--word_size', type=int, help="kmer_size")

    parser.add_argument('-m', '--mode', type=str, choices=['greedy', 'fast', 'deep'], 
    					default='fast', help="The algorithm to generate the kmer strings for hashing")

    return parser.parse_args(cmdline_params)


class Spire(object):
	
	"""This is where the magic happens"""

	CALCULATORS = [SingleSketchCalculator, ExhaustivePairSketchCalculator]

	def __init__(self, word_size, n_hash_functions, mode):
		
		self._word_size = word_size
		self._n_hash_functions = n_hash_functions
		self._calculator = self.setMode(mode)
		self.ref_signature = None

	def setMode(self, mode):
		for calculator in self.CALCULATORS:
			if calculator.canCalculate(mode):
				return calculator(self._word_size, self._n_hash_functions, MinHash)

	def loadRefSignature(self, ref):
		self.ref_signature = self._calculator.run(ref)

	def getJaccard(self, query_sequence):
		query_signature = self._calculator.run(query_sequence)
		return self.approximateJaccard(self.ref_signature, query_signature)

	def approximateJaccard(self, sig1, sig2):
		"""
		Calculate an approximate Jaccard index based on minhash Sketchs
		TODO: Assert that sig1 and sig2 were calculated with the same hash functions, otherwise
		the approimation is invalid.
		:param list(int) sig1:
		:param sig2
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
	print jaccard
