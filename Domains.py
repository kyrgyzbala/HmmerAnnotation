'''
Created on Jun 2, 2014

@author: sanjarbek
'''

import params as p

class Coordinate(object):
    def __init__(self,pFrom=None,pTo=None):
        self.pFrom = pFrom
        self.pTo = pTo
    
    def __str__(self):
        return "Coordinate: (from,to)=(%d, %d)"%(self.pFrom, self.pTo)

class DomainSceleton(object):
    
    def extend(self,seq):
        self.ext_coordinate = Coordinate()
        total_seq_length = len(seq.seq)
        right = 0
        stop = False
        while not stop and right<p.max_extend_range and self.coordinate.pTo+right < total_seq_length:
            if seq.seq[self.coordinate.pTo-1+right]=='*':
                stop = True
                self.ext_coordinate.pTo = self.coordinate.pTo+right
            else:
                right+=1
            if self.coordinate.pTo-1+right==len(seq.seq)-1: # if it reaches the end
                stop = True
        if not self.ext_coordinate.pTo:
            self.ext_coordinate.pTo = self.coordinate.pTo
        
        self.ext_coordinate.pFrom = None
        inc = 0
        stop = False
        leftmost_M,leftmost_V,previous_stop = None,None,None
        #search for stop codon and M (or V) to the left
        #At first it seems that meeting M might be sufficient reason to stop. But it's a trap!
        #M can be in the middle of a protein as well.
        
        while not stop:
            if seq.seq[self.coordinate.pFrom-1-inc] == 'M':
                leftmost_M = self.coordinate.pFrom-inc
            elif seq.seq[self.coordinate.pFrom-1-inc] == 'V':
                leftmost_V = self.coordinate.pFrom-inc
            elif seq.seq[self.coordinate.pFrom-1-inc]=='*':
                stop = True
                previous_stop = self.coordinate.pFrom-inc
            elif self.coordinate.pFrom-inc<=0:
                stop = True
            inc+=1
        
        if leftmost_M:
            left_edge = leftmost_M
        elif leftmost_V:
            left_edge = leftmost_M
        elif previous_stop:
            left_edge = previous_stop
        else:
            left_edge = self.coordinate.pFrom
        
        #search for M to the right (this is for the case when hit overspans the actual sequence's start point)
        inc = 1
        stop = False
        right_edge = None
        while not stop and inc<10:
            if seq.seq[self.coordinate.pFrom+inc]=='M':
                right_edge = self.coordinate.pFrom+inc
                stop = True
            inc+=1
        #decide which one to set as border
        if right_edge and left_edge:
            if abs(left_edge-self.coordinate.pFrom)/abs(right_edge-self.coordinate.pFrom)>p.allowed_ratio_for_edges:
                self.ext_coordinate.pFrom = right_edge
            else:
                self.ext_coordinate.pFrom = left_edge
        elif right_edge and not left_edge:
            self.ext_coordinate.pFrom = right_edge
        elif not right_edge and left_edge:
            self.ext_coordinate.pFrom = left_edge
        elif not right_edge and not left_edge:
            self.ext_coordinate.pFrom = self.coordinate.pFrom
        
        self.sequence = seq.seq[self.coordinate.pFrom - 1:self.coordinate.pTo-1]
        self.e_sequence = seq.seq[self.ext_coordinate.pFrom-1:self.ext_coordinate.pTo-1]

    def inc_coordinates(self, inc):
        """Written for turning coordinate within a tile into global protein coordinate.
        inc is tile.pFrom
        Not to be confused with scale_back"""
        self.coordinate.pFrom += inc
        self.coordinate.pTo   += inc
        self.ext_coordinate.pFrom += inc
        self.ext_coordinate.pTo   += inc

    def scale_back(self,length):
        """Convert protein coordinates to nucleotide coordinates pr_cr*3+frame.
        Be sure to call it while you don't lose the tile coordinates (add tile.pFrom to all the coordinates using inc_coordinates), 
        in case if you are running on tiles."""
        if self.frame>0:
            self.sc_coordinate = Coordinate(self.coordinate.pFrom*3+self.frame,self.coordinate.pTo*3+self.frame)
            self.sc_ext_coordinate = Coordinate(self.ext_coordinate.pFrom*3+self.frame,self.ext_coordinate.pTo*3+self.frame)
        else:
            self.sc_coordinate = Coordinate(length-self.coordinate.pTo*3+1,length-self.coordinate.pFrom*3+1)
            self.sc_ext_coordinate = Coordinate(length-self.ext_coordinate.pTo*3+1,length-self.ext_coordinate.pFrom*3+1)
    
    def __lt__(self,other):
        if self.ext_coordinate:
            return self.sc_ext_coordinate.pFrom < other.sc_ext_coordinate.pFrom
        else:
            return self.sc_coordinate.pFrom < other.sc_coordinate.pFrom

    def __str__(self):
        return self.name

class HitRegion(DomainSceleton):
    # or envelope from hmmer domtblout
    def __init__(self,hmmer_line):
        parts = hmmer_line.split()
        self.frame = int(parts[-1])
        pFrom, pTo = int(parts[p.ecc[0]-1]), int(parts[p.ecc[1]-1])
        
        self.coordinate = Coordinate(pFrom,pTo)        
        # ext_coordinate = None denotes that domain extension did not happen yet.        
        self.ext_coordinate = None
        # sc_coordinate and sc_ext_coordinate means the scale back operation is not carried out yet. 
        self.sc_coordinate = None
        self.sc_ext_coordinate = None
        
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
