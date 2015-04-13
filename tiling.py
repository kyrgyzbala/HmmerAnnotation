import tools
from TileClasses import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import os
import params as p


def decompose(fna_file,tiling_length,tiling_overlap):
	fna_seq = SeqIO.read(fna_file,'fasta')
	total_length = len(fna_seq.seq)
	pivot = 0
	tiles = []

	while pivot<total_length:
		if pivot+tiling_length<total_length:
			curSeq = Seq(str(fna_seq.seq[pivot:pivot+tiling_length]))
			curSeqRecord = SeqRecord(curSeq,id='',description='')
			curCoordinate = Coordinate(pivot,pivot+tiling_length)
		else:
			curSeq = Seq(str(fna_seq.seq[pivot:]))
			curSeqRecord = SeqRecord(curSeq,id='',description='')
			curCoordinate = Coordinate(pivot,total_length-1)
		tiles.append(Tile(curSeqRecord,curCoordinate))
		pivot += tiling_length - tiling_overlap
	return tiles

def search_tile(tile,hmmProfile,evalue_thr):
	
	seq_file = 'tile.fna'
	SeqIO.write(tile.seq,seq_file,'fasta')
	regions = tools.get_hit_regions(seq_file,hmmProfile,evalue_thr)
	os.remove(seq_file)
	return regions

def generate_annotations(fna_file,hmmProfile,tiling_length=p.tiling_length,tiling_overlap=p.tiling_overlap,evalue_thr=p.dom_evalue_thr):
	tiles = decompose(fna_file,tiling_length,tiling_overlap)
	out_lines = []
	all_regions=[]
	source = os.path.basename(fna_file).split('.')[0]
	name = hmmProfile.name if hmmProfile.name else hmmProfile.hmm_file_path.split('.')[0]
	for tile in tiles:
		result = search_tile(tile,hmmProfile,evalue_thr)
		out_fmt = "%s\t%d\t%s\t%s\t%s\t%s\t%s"
		tmp_id = 1
		if result.region:
			start = result.region[0]+tile.coordinate.start
			end = result.region[1]+tile.coordinate.start
			strand = result.region[2]
			out_str = out_fmt%(source,tmp_id,strand,start,end,name,'hmmer hit from %s'%name)
			tmp_id+=1
			out_lines.append(out_str)
		elif result.regions:
			for i in range(len(result.regions)):
				result.regions[i].start+=tile.coordinate.start
				result.regions[i].end+=tile.coordinate.start
			all_regions+=result.regions

	# for region in all_regions:
	# 	print region.start, region.end, region.strand

	mergedRegions = result.aggregate(all_regions)

	for region in mergedRegions:
		start = region.start
		end = region.end
		strand = region.strand
		out_str = out_fmt%(source,tmp_id,strand,start,end,name,'hmmer hit from %s'%name)
		tmp_id+=1
		out_lines.append(out_str)

	return out_lines