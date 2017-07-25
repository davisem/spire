__author__ = "Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""

import unittest
from sketch_calculators import SingleSketchCalculator, ExhaustivePairSketchCalculator

class SketchCalculatorTests(unittest.TestCase):

	def test_getSketch(self):
		"""Test that we can get a sketch"""
		test_seq = "AAATTACCCGGG"
		expect_length = 100

		Calculator = SingleSketchCalculator(3, expect_length)
		sketch = Calculator.getSketch(test_seq)
		self.assertEquals(len(sketch), expect_length)

	def test_setWordSize(self):
		"""Test that the word size for a pair calculator is double the word size of a single calculator"""
		pair_calculator = ExhaustivePairSketchCalculator(2, 1)
		single_calculator = SingleSketchCalculator(2, 1)
		self.assertEqual(pair_calculator.word_size * 2, single_calculator.word_size)
