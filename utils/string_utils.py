#!/usr/bin/env python

__author__ = "Pattrick Hagarty, Eric Davis"
__copyright__ = ""
__credits__ = ["Eric Davis"]
__version__ = ""
__maintainer__ = "Eric Davis"
__email__ = "emdavis48@gmail.com"
__status__ = ""
__modname__ = "string_utils.py"

import numpy as np


def JaccardViaSets(words1, words2):
	"""Calculates a Jaccard index using set logic"""

	return float(len(set(words1) & set(words2))) / len(set(words1) | set(words2))


def CreateWindow(wl):
	"""Creates a random DNA string of a certain size"""
	window = ""
	r = np.random.randint(4, size=wl)
	for i in range(wl):
		if r[i] == 0:
			window += 'G'
		elif r[i] == 1:
			window += 'T'
		elif r[i] == 2:
			window += 'C'
		else:
			window += 'A' 
	return window

def CreateKMERPairs ( w, k ):
    kmerpairs = {}
    for i in range ( len(w) - k):
        for j in range (i+1,len(w) - k ):
            a = w[i:i+k]
            b = w[j:j+k]
            d = j - i
            ab= a+","+b+":"
            if ab in kmerpairs:
                kmerpairs[ab] = min(d,kmerpairs[ab])
            else :
                kmerpairs[ab] = d
    return kmerpairs
    
def PerturbWindow(original_window, prob):
	"""Perturbs a given window given some probability"""

	window = ""
	r = np.random.rand(len(original_window))
	for i in range ( len(original_window)):
		if r[i] < prob :
			rr = np.random.randint(6)
			if rr == 0:
				window += 'G'
			elif rr == 1:
				window += 'T'
			elif rr == 2:
				window += 'C'
			elif rr == 3:
				window += 'A'
			elif rr == 4:
				rrr = np.random.randint(4)
				if rrr == 0:
					window += 'G'
				elif rrr == 1:
					window += 'T'
				elif rrr == 2:
					window += 'C'
				elif rrr == 3:
					window += 'A'
				i -= 1
		else:
			window += original_window[i]
	return window
