import params as p
# from DomainClasses import *
from Domains import HitRegion
from Bio import SeqIO
from translate import translate
import subprocess as sp
import os
from HmmProfileClasses import HmmResult
import sys

def get_hit_coordinates(file_name,trans_seqs,length,evalue_thr=p.dom_evalue_thr):
	hmm_hits = [l for l in open(file_name).readlines() if not l.startswith('#') and l.strip()!='']
	hmm_hits = [l for l in hmm_hits if float(l.split()[p.dec-1])<=float(evalue_thr)]
	domains = [SingleDomain(l) for l in hmm_hits]
	[d.extend(trans_seqs[p.frames.index(d.frame)],length) for d in domains]
	# for d in domains:
	# 	print d.start,d.end,d.sequence
	# 	print d.e_start,d.e_end,d.e_sequence
	domains.sort()
	return domains

def parse_domains(file_name,trans_seqs,evalue_thr=p.dom_evalue_thr):
	
	hmm_hits = [l for l in open(file_name).readlines() if not l.startswith('#') and l.strip()!='']
	hmm_hits = [l for l in hmm_hits if float(l.split()[p.dec-1])<=float(evalue_thr)]
	domains = [HitRegion(l) for l in hmm_hits]
	# Extension procedure can be taken out of this procedure as well.
	# I put it here for reducing IO by reading (Bio.SeqIO) the trans_seq once, and using it here.
	[d.extend(trans_seqs[p.frames.index(d.frame)]) for d in domains]
	
	return domains

def get_hit_regions_new(seq_file,hmmProfile):
	
	seq = SeqIO.read(seq_file,'fasta')
	global_length = len(seq.seq.tostring())
	seqs = translate(seq)
	hmm_out_file = '%s.out'%seq_file.split('.')[0]
	all_translated = '%s_6f.faa'%seq_file.split('.')[0]
	SeqIO.write(seqs,open(all_translated,'w'),'fasta')
	sp.call(['hmmsearch','--domtblout',hmm_out_file,hmmProfile.file_path,all_translated],stdout=sp.PIPE)
	
	hits = parse_domains(hmm_out_file, seqs, global_length)
	
	os.remove(hmm_out_file)
	os.remove(all_translated)
	return hits

def get_hit_regions(seq_file,hmm_profile,evalue_thr):
	
	seq = list(SeqIO.parse(seq_file,'fasta'))[0]
	global_length = len(seq.seq.tostring())
	seqs = translate(seq)
	hmm_out_file = '%s.out'%seq_file.split('.')[0]
	all_translated = '%s_6f.faa'%seq_file.split('.')[0]
	SeqIO.write(seqs,open(all_translated,'w'),'fasta')
	sp.call(['hmmsearch','--domtblout',hmm_out_file,hmm_profile.file_path,all_translated],stdout=sp.PIPE)
	hmmResult = HmmResult(hmm_profile,hmm_out_file)
	hits = get_hit_coordinates(hmm_out_file, seqs, global_length,evalue_thr)
	if hits:
		hmmResult.hits = hits
		hmmResult.validate()
	else:
		hmmResult.regions = []
	os.remove(hmm_out_file)
	os.remove(all_translated)
	return hmmResult

def load_ptt_dict(ptt_coordinate_file):
	crd_dict = {}
	for l in open(ptt_coordinate_file).readlines():
		terms = l.split(',')
		if terms[0] in crd_dict:
			crd_dict[terms[0]].append((terms[1],int(terms[2]),int(terms[3])))
		else:
			crd_dict[terms[0]] = [(terms[1],int(terms[2]),int(terms[3]))]
	return crd_dict

def generate_annotations(sequence_file, hmm_profile, evalue_thr):
	
	annotations = get_hit_regions(sequence_file, hmm_profile, evalue_thr)
	source_name = os.path.splitext(os.path.basename(sequence_file))[0]
	
	annotation_lines = []
	tmp_id = 1
	str_fmt = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t"
	for r in annotations.regions:
		l = str_fmt%(source_name,tmp_id,r.strand,r.start,r.end,'-',hmm_profile.source)
		annotation_lines.append(l)
		tmp_id+=1
	return annotation_lines

def sign(orf):
	if type(orf)=='str':
		orf = int(orf)
	return '+' if orf>0 else '-'

def get_real_cords(src,domain,ptt_dict):
	res_crd = ()
	f_p,t_p,str_p = domain.scaled_coords[0],domain.scaled_coords[1],sign(domain.frame)
	dist = 1e20 #some very big number
	for tmp_crd in ptt_dict[src]:
		str_a,f_a,t_a = tmp_crd[0],tmp_crd[1],tmp_crd[2]
		if str_a==str_p:
			tmp_dist = abs(f_p-f_a)+abs(t_p-t_a)
			if tmp_dist<dist:
				res_crd = (f_a,t_a)
				dist = tmp_dist
	return res_crd


def overlapping(first,second):
	#Arguments 'first' and 'second' must be instances of Coordinate class	
	dist1 = first.pFrom-second.pTo
	dist2 = second.pFrom-first.pTo
	if dist1*dist2>=0:
		return True
	return False
