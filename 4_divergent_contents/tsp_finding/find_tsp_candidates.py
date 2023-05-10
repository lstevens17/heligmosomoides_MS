from statistics import median
import ete3
import sys

# usage:
#       python3 find_tsps.py <treefile> <boostrapthreshold>

# set list of haplotypes which we'll check for monophyly
hb_haps = ['nxHelBake1.alternate', 'nxHelBake2.primary', 'nxHelBake3.alternate', 'nxHelBake3.primary', 'nxHelBake2.alternate', 'nxHelBake1.primary']
hp_haps = ['ngHelPoly1.alternate', 'ngHelPoly2.primary', 'ngHelPoly2.alternate', 'ngHelPoly1.primary']

# set outgroup 
outgroup = "h_mixtum"

# read in treefile into dict and root
with open(sys.argv[1], 'r') as treefile:
    tree_dict = {}
    for line in treefile:
        id, nwk = line.rstrip("\n").split("\t")[0], line.rstrip("\n").split("\t")[1]
        tree = ete3.Tree(nwk)
        tree.set_outgroup(outgroup)
        tree_dict[id] = tree

# set bootstrap threshold
bs_threshold = int(sys.argv[2])

for id, tree in tree_dict.items(): 
    # first identify suspicious leaves
    leaves = tree.get_leaf_names() # get a list of leaf names
    leaves.remove(outgroup) # remove the outgroup
    dist_list = []
    sus_list = []
    for leaf in leaves: # for every leaf
        dist_list.append(tree.get_distance(outgroup, leaf)) # store its distance to outgroup
    for leaf in leaves: # for every leaf
        if tree.get_distance(outgroup, leaf) > median(dist_list) * 2: # if its distance to root is 2x median
            sus_list.append(leaf) # add it to the supicious list
    prune_list = list(set(leaves) - set(sus_list))  # get a list of non-suspicious leafs
    prune_list.append(outgroup) # add outgroup to it
    tree.prune(prune_list) # prune the tree to include only non-suspicious sequences
    tsp_candidate = "-"
    hb1_hb1_status = "-"
    # now identity trees containing mixed clades
    for node in tree.traverse("preorder"): # traverse tree
        if node.is_root() == False and node.is_leaf() == False: # ignore root and leaf nodes
            leaves = node.get_leaf_names()
            prefix_list = []
            if len(leaves) != len(prune_list) - 1: # if we're not at the LCA of all hb/hp haps
                for leaf in leaves: # for every leaf
                    prefix = leaf[0:9]
                    prefix_list.append(prefix) # store the prefix (e.g. nxHelBake) in a list
                if len(set(prefix_list)) == 2: # if both nxHelBake and ngHelPoly in list
                    if node.support >= bs_threshold: # and clade has support above bootstrap threshold
                        tsp_candidate = "TSP" # record as TSP candidate
                        if "nxHelBake1.primary" in leaves and "ngHelPoly1.primary" in leaves: # if the clade has both HB1 and HP1
                            hb1_hb1_status = "Y" # record this
    if len(sus_list) == 0:
        sus_list.append("-")
    # print output file (tree id, TSP status, suspicious leaves, pruned tree)
    print("%s\t%s\t%s\t%s\t%s" % (id, tsp_candidate, hb1_hb1_status,  ",".join(sus_list), tree.write()))

    
