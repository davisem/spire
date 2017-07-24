TODO: implement this

class Winnower(object):
	def __init__(self):
		pass


	def storeSketch(self, seq, Sketch):
		for i, min_mer in enumerate(Sketch):
			self._Sketch_cache[i][min_mer].append(seq)

	def queryRead(self, Sketch):
		for i, min_mer in enumerate(Sketch):
			yield self._Sketch_cache[i][min_mer]

	def getSimilarReads(self, read):
		
		similar_reads = []
		for q_read_list in self.queryRead(read.sketch):
			for q_read in q_read_list:
				if q_read != read.seq:
					similar_reads.append(q_read)
		
		return Counter(similar_reads)