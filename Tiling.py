from Bio import SeqIO
from BioClasses import Gene
import os
import tools
import params as p
import sys

class Tile(Gene):
	def __init__(self, source, pFrom, pTo):
		super(Tile, self).__init__(source, pFrom, pTo)
	
	def __str__(self):
		return  "Tile of %s,Coords: %d:%d"%(self.src, self.pFrom, self.pTo)



def decompose(fna_file,tiling_length,tiling_overlap):
	fna_seq = SeqIO.read(fna_file,'fasta')
	src_name = os.path.basename(fna_file).split('.')[0]
	total_length = len(fna_seq.seq)
	pivot = 1
	tiles = []
	
	while pivot<total_length:
		if pivot+tiling_length<total_length:
			curTile = Tile(source=src_name, pFrom=pivot, pTo=pivot+tiling_length)
		else:
			curTile = Tile(source=src_name, pFrom=pivot, pTo=total_length-1)
		
		curTile.set_sequence(fna_seq)
		tiles.append(curTile)
		pivot += tiling_length - tiling_overlap
	
	return tiles, total_length

def search_tile(tile,hmmProfile,evalue_thr):
	
	seq_file = 'tile.fna'
	SeqIO.write(tile.seq,seq_file,'fasta')
	hits = tools.get_hit_regions_new(seq_file,hmmProfile)
	os.remove(seq_file)
	return hits

def get_hit_regions(fna_file,hmmProfile,tiling_length=p.tiling_length,tiling_overlap=p.tiling_overlap,evalue_thr=p.dom_evalue_thr):
	
	tiles, global_seq_len = decompose(fna_file,tiling_length,tiling_overlap)
	
	all_hits = []
	for tile in tiles:
		hits = search_tile(tile,hmmProfile,evalue_thr)
		if hits:
			[hit.inc_coordinates(tile.coordinate.pFrom) for hit in hits]
			all_hits+=hits
	
	[hit.scale_back(global_seq_len) for hit in all_hits]
	all_hits.sort()
	
	# I expect there should be no overlapping hits after extending them.
	prev_hit = None
	for hit in hits:
		if not prev_hit:
			prev_hit=hit
			continue
		are_overlapping = True if tools.overlapping(prev_hit.sc_ext_coordinate, hit.sc_ext_coordinate) else False
		assert not (are_overlapping and prev_hit.frame==hit.frame), "Overlapping hit found."
	
	return all_hits
	
	
	

