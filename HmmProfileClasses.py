import params as p
from DomainClasses import *
import DomainClasses as dc
import os

def in_proximity(first, second):
	(first_start,first_end) = first.scaled_coords
	(second_start,second_end) = second.scaled_coords
	tst_diff = second_start - first_end
	dist = tst_diff if tst_diff>0 else first_start - second_end
	if dist<10000:
		return True
	else:
		return False

def overlapping(first,second):
	try:
		(first_start,first_end) = first.scaled_coords
		(second_start,second_end) = second.scaled_coords
	except:
		(first_start,first_end) = first.start,first.end
		(second_start,second_end) = second.start,second.end
	dist1 = first_start - second_end
	dist2 = second_start - first_end
	if dist1*dist2>=0:
		return True
	return False

class HmmProfile:
	def __init__(self, hmm_file_path, profile_type=None, order=[], name='not-named'):
		self.file_path = hmm_file_path
		self.source = os.path.splitext(os.path.basename(hmm_file_path))[0]
		self.raw_lines = open(hmm_file_path).readlines()
		self.special_case=profile_type
		self.order=order
		self.name=name

class HmmResult:
	def __init__(self, HmmProfile=None, hmm_out_file_path=None):
		self.HmmProfile = HmmProfile
		self.file_path = hmm_out_file_path
		self.hits = []
		self.region = None
	
	def validate(self):

		if self.HmmProfile.special_case=='LuxR':
			names = [d.name for d in self.hits]
			if sorted(names)==p.luxr_domains:
				region = dc.aggregate(self.hits)
			else:
				region = None
			self.region = region
		else:
			self.regions = self.aggregate(self.hits)

	def aggregate(self,region_array):
		if len(region_array)==0:
			return []
		groups = []
		tmp_group = []
		prev_dom = None
		for i in range(len(region_array)):
			domain = region_array[i]
			if not prev_dom:
				prev_dom = domain
				tmp_group.append(domain)
				if len(region_array)==1:
					groups.append(tmp_group)
				continue
			try:
				prev_frame = prev_dom.frame
				cur_frame = domain.frame
			except:
				prev_frame = 1 if prev_dom.strand=='+' else -1
				cur_frame = 1 if domain.strand=='+' else -1

			if overlapping(prev_dom,domain) and prev_frame*cur_frame>0:
				tmp_group.append(domain)
				prev_dom = domain
				if i == len(region_array)-1:
					groups.append(tmp_group)
					continue
			else:
				groups.append(tmp_group)
				tmp_group = [domain]
				prev_dom = domain
			if i == len(region_array)-1:
					groups.append(tmp_group)
		# for group in groups:
		# 	print 'The group:'
		# 	for d in group:
		# 		print '  ', d.name,d.start, d.end, d.e_start, d.e_end, d.scaled_coords, d.evalue

		regions = [AggregatedRegion(group,self.HmmProfile.special_case) for group in groups]
		regions = [r for r in regions if r.valid]
		return regions