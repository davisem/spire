__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"


import cProfile
import unittest

from pstats import Stats

from spire import Spire
from utils.string_utils import CreateWindow, CreateKMERPairs, PerturbWindow

class SpireTests(unittest.TestCase):
	
	def setUp(self):
		"""Enable profiling"""
		
		self.pr = cProfile.Profile()
		self.pr.enable()
		print "\n<<<---"

	def tearDown(self):
		"""Report profiling results"""
		
		p = Stats(self.pr)
		p.strip_dirs()
		p.sort_stats('cumtime')
		p.print_stats()
		print "\n--->>>"

	def test_approximateJaccardGood(self):
		"""Test that the approximate jaccard is 1 at max identity"""

		sig1 = [20, 20, 10]
		sig2 = [20, 20, 10]
		spire = Spire(2, 3, "fast")
		jaccard = spire.approximateJaccard(sig1, sig2)
		
		self.assertEquals(1, jaccard)

	def test_approximateJaccardBad(self):
		"""Test that the approximate jaccard is 0 at max disparity"""
		
		sig1 = [100, 200, 30]
		sig2 = [20, 20, 10]
		spire = Spire(2, 3, "fast")
		jaccard = spire.approximateJaccard(sig1, sig2)
		
		self.assertEquals(0, jaccard)

	def test_jaccard_single(self):
		"""Test that single mode returns a sane jaccard value between 0 and 1"""
		
		test_seq = CreateWindow(1000)
		random_seq = CreateWindow(1000)
		q_test_seq = PerturbWindow(test_seq, 0.3)
		spire = Spire(5, 400, "fast")
		
		spire.loadRefSignature(test_seq)
		print(spire.getJaccard(q_test_seq))
		print("random seq {}".format(spire.getJaccard(random_seq)))
		self.assertTrue(0 < spire.getJaccard(q_test_seq) < 1)

	def test_jaccard_deep(self):
		"""Test that deep mode returns a sane jaccard value between 0 and 1"""
		
		test_seq = CreateWindow(1000)
		random_seq = CreateWindow(1000)
		q_test_seq = PerturbWindow(test_seq, 0.3)
		spire = Spire(5, 400, "deep")
		
		spire.loadRefSignature(test_seq)
		print(spire.getJaccard(q_test_seq))
		print("random seq {}".format(spire.getJaccard(random_seq)))
		self.assertTrue(0 < spire.getJaccard(q_test_seq) < 1)

