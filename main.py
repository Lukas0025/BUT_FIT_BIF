##
# @author Lukáš Plevač <xpleva07@vutbr.cz>
# @date 12.5.2023
# Semestral project to BIF on BUT FIT
# Create ancestrals amino acides seqvences for Phylogenetic tree
#

from Bio import Phylo
from Bio import AlignIO
import pandas as pd

DEBUG           = False
TREE_FILE       = "tree.tre"
MSA_FILE        = "msa.fasta"
ANCESTRALS_FILE = "ancestrals.csv"
OUT_DIR         = "output"

##
## FILE LOAD PART
##

# first load phylo genetic tree
PhyloTree = Phylo.read(TREE_FILE, "newick")

# second load fasta
FastaAlign = AlignIO.read(MSA_FILE, "fasta")

# third load ancestrals
ancestrals = pd.read_csv(ANCESTRALS_FILE, sep = ",")


# create dict for simple fasta indexing by fasta ID
FastaAlignDict = {}
for fasta in FastaAlign:
    FastaAlignDict[fasta.id] = fasta

##
## Check if tree file coresponding with ancestrals
##

## get all nodes from ancestrals
a_nodes = []
for node in ancestrals['node']:
    if node not in a_nodes:
        a_nodes.append(node)

a_nodes.sort()

## get all nodes from tree
t_nodes = []
for clade in PhyloTree.find_clades():
    if clade.confidence is not None:
        t_nodes.append(clade.confidence)

t_nodes.sort()

if a_nodes != t_nodes:
    print("Error not same nodes in tree and ancestrals files")
    print("Nodes in tree file:       {}".format(t_nodes))
    print("Nodes in ancestrals file: {}".format(a_nodes))
    exit(1)

##
## Create acentral seqvence without spaces
## using posterior probability form csv file
##

seq = {}

for node in a_nodes:
    seq[node] = {}
    
    ac_node   = ancestrals[ancestrals['node'] == node]
    ac_node.reset_index()

    for _, row in ac_node.iterrows():
        amino = row[['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']]
        amino = pd.to_numeric(amino, errors ='coerce').fillna(0).astype('float')
        
        seq[node][row['position'] - 1] = amino.idxmax(axis=0)

##
## Space placement part
##

##
# Get all leaves nodes for clade and distance to it
#
# @param clade clade to find all leafs nodes
# @return array of dict of leaf clade and total_branch_len from param clade
def getLeaves(clade):
    opend  = [{'clade': clade, 'total_len': 0}]
    closed = []

    while len(opend) > 0:
        parrent = opend.pop()
        if parrent['clade'].name is not None:
            closed.append(parrent)
        else:
            for kid in parrent['clade']:
                opend.append({'clade': kid, 'total_len': parrent['total_len'] + kid.branch_length})

    return closed

## Naive space place by leaves but not length
## commented not used
'''
for clade in PhyloTree.find_clades():
    if clade.confidence is not None:
        leaves      = getLeaves(clade)
        leaves_seqs = [FastaAlignDict[leaf['clade'].name] for leaf in leaves]
        
        for pos in seq[clade.confidence].keys():
            spaces = 0
            for leave_seq in leaves_seqs:
                if leave_seq[pos] == "-":
                    spaces += 1

            space_prob = spaces / len(leave_seq)

            if (space_prob >= 0.5):
                seq[clade.confidence][pos] = '-'
'''

## space placement by leaves and count with length of branches
for clade in PhyloTree.find_clades():
    if clade.confidence is not None:
        leaves      = getLeaves(clade)
        leaves_seqs = [{'seq': FastaAlignDict[leaf['clade'].name], 'len': leaf['total_len']} for leaf in leaves]
        
        for pos in seq[clade.confidence].keys():
            space_prob  = 0
            len_sum     = 0

            for leave_seq in leaves_seqs:
                len_sum += leave_seq['len']

            for leave_seq in leaves_seqs:
                if leave_seq['seq'][pos] == "-":
                    space_prob += leave_seq['len'] / len_sum

            if (space_prob >= 0.5):
                seq[clade.confidence][pos] = '-'


##
## Save seq to files
##

# first convert seq dict to string
for name, sq in seq.items():
    str_seq = ['-'] * len(sq)
    
    for index, amino in sq.items():
        str_seq[index] = amino

    seq[name] = "".join(str_seq)

    if DEBUG:
        print("{}: {}".format(name, seq[name]))

# save to files
for name in seq.keys():
    text_file = open("{}/node_{}.fas".format(OUT_DIR, name), "wt")
    text_file.write(seq[name])
    text_file.close()
