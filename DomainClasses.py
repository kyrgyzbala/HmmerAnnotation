import params as p


class DomainSceleton:
	
	def extend(self,seq,length):
		right = 0
		stop = False
		while not stop and right<100 and self.end+right < len(str(seq.seq)):
			if seq.seq[self.end-1+right]=='*':
				stop = True
				self.e_end = self.end+right
			else:
				right+=1
			if self.end-1+right==len(seq.seq)-1:
				stop = True
		if not self.e_end:
			self.e_end = self.end
		self.e_start = None
		inc = 0
		stop = False
		leftmost_M,leftmost_V,previous_stop = None,None,None
		#search for stop codon and M (or V) to the left
		while not stop:
			if seq.seq[self.start-1-inc] == 'M':
				leftmost_M = self.start-inc
			elif seq.seq[self.start-1-inc] == 'V':
				leftmost_V = self.start-inc
			elif seq.seq[self.start-1-inc]=='*':
				stop = True
				previous_stop = self.start-inc
			elif self.start-inc<=0:
				stop = True
			inc+=1
		if leftmost_M:
			left_edge = leftmost_M
		elif leftmost_V:
			left_edge = leftmost_M
		elif previous_stop:
			left_edge = previous_stop
		else:
			left_edge = self.start
		#search for M to the right (this is for the case when hit overspens the actual sequence)
		inc = 1
		stop = False
		right_edge = None
		while not stop and inc<10:
			if seq.seq[self.start+inc]=='M':
				right_edge = self.start+inc
				stop = True
			inc+=1
		#decide which one to set as border
		if right_edge and left_edge:
			if abs(left_edge-self.start)/abs(right_edge-self.start)>100:
				self.e_start = right_edge
			else:
				self.e_start = left_edge
		elif right_edge and not left_edge:
			self.e_start = right_edge
		elif not right_edge and left_edge:
			self.e_start = left_edge
		elif not right_edge and not left_edge:
			self.e_start = self.start		
		self.scale_back(length)
		self.sequence   = seq.seq[self.start - 1:self.end-1]
		self.e_sequence = seq.seq[self.e_start-1:self.e_end-1]

	def scale_back(self,length):
		if self.frame>0:
			self.scaled_coords = (self.start*3+self.frame,self.end*3+self.frame)
			self.e_scaled_coords = (self.e_start*3+self.frame,self.e_end*3+self.frame)
		else:
			self.scaled_coords = (length-self.end*3+1,length-self.start*3+1)
			self.e_scaled_coords = (length-self.e_end*3+1,length-self.e_start*3+1)

	def __lt__(self,other):
		if self.e_scaled_coords:
			return self.e_scaled_coords[0] < other.e_scaled_coords[0]
		else:
			return self.scaled_coords[0] < other.scaled_coords[0]

	def __str__(self):
		return self.name

class SingleDomain(DomainSceleton):
	
	def __init__(self,hmmer_line):
		parts = hmmer_line.split()
		self.frame = int(parts[-1])
		self.start = int(parts[p.ecc[0]-1])
		self.end = int(parts[p.ecc[1]-1])
		self.e_start = None
		self.e_end = None
		self.name = parts[3]
		self.evalue = parts[p.dec-1]

class SuperDomain(DomainSceleton):
	def __init__(self, SingleDomains):
		self.SingleDomains = SingleDomains.sort()

class AggregatedRegion:
	def __init__(self, regions,profile_type):
		self.regions = regions
		self.special_case = profile_type
		self._aggregate()
		self.sequences = []
	def _aggregate(self):
		if self.special_case=='general':
			self.start, self.end, self.strand, self.valid, self.sequences = aggregate(self.regions)
		elif self.special_case=='LuxR':
			self.start, self.end, self.strand, self.valid, self.sequences = aggregate_luxr(self.regions)
	def __str__(self):
		if self.special_case:
			return 'AggregatedRegion '+self.special_case
		else:
			return 'AggregatedRegion general'


def aggregate_luxr(regions):
	domains = set(r.name for r in regions)
	if len(domains)>=2 and 'GerE' in domains:
		return aggregate(regions)
	else:
		return 0,0,0,False,[]

def aggregate(regions):
	try:
		start = min(r.e_scaled_coords[0] for r in regions)
		end = max(r.e_scaled_coords[1] for r in regions)
		strand = '+' if regions[0].frame>0 else '-'
		sequences = [r.e_sequence.tostring() for r in regions]
	except:
		start = min(r.start for r in regions)
		end = max(r.end for r in regions)
		strand = regions[0].strand
		sequences = ''
	
	return start,end,strand,True,sequences
