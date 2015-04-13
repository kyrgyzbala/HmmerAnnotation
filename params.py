
frames = [1,2,3,-1,-2,-3]
ecc = (20,21) # envelope coordinate column numbers from hmmer tabular output.
dec = 12 #column which contains evalues for domain in the domain mode tabular output
dom_evalue_thr = 1e-10
domain_dist = 1000 #distance threshold for merging domains

max_extend_range = 100

# If the hit's starting point overspans the starting point of actual domain, 
# the starting M(or V) will be located on the right side. Since the straightforward 
# logic is to search for starting point to the left, it will make large margin mistake.
# Therefore, there's an option of searching towards the right. And if there's an M on the right side 
# very close (defined by max_M_inner_search_range) to the starting point of hit, then most likely it 
# is the actual starting point and in that case you shrink the domain, rather than extending.
# The decision is made in terms of ratio of distances of the possible extensions coordinates from the left
# side over right side. if left_stop/right_stop > allowed_ratio_for_edges, then it chooses the right side. 
max_M_inner_search_range = 10
allowed_ratio_for_edges = 100

tiling_length = 5000
tiling_overlap = 1000

luxr_domains = ['IPR000792','IPR005143']