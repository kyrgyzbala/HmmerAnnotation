'''
Created on Jun 1, 2014

@author: sanjarbek
'''

import Tiling
import Profile
import sys


if __name__=='__main__':
    fna_file = '/home/sanjarbek/Data/NCBI/complete/Burkholderia_thailandensis_E264_uid58081/NC_007650.fna'
    profile_file = '../files/profiles/IPRall_R.hmm'
    hmmProfile = Profile.HmmProfile(profile_file)
    
    hits = Tiling.get_hit_regions(fna_file, hmmProfile)
    outfmt="%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
    print "\t".join(["Pname","Frm","From","To","sc_From","sc_To","e_From","e_To","e_sc_From","e_sc_To"])
    for hit in hits:
        
        print outfmt%(hit.name,hit.frame, hit.coordinate.pFrom, hit.coordinate.pTo, \
                      hit.sc_coordinate.pFrom, hit.sc_coordinate.pTo, \
                      hit.ext_coordinate.pFrom, hit.ext_coordinate.pTo, \
                      hit.sc_ext_coordinate.pFrom, hit.sc_ext_coordinate.pTo )
    
    