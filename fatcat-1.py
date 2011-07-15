#!/usr/bin/env python

import os
import ete2a1
from ete2a1 import PhyloTree
from Bio import SeqIO
from optparse import OptionParser
import sys

def insert_sequence(pt,highest_node,seq, base,aln_name):
    print 'insert_sequence: Before'
    print pt
    seq_handle = open(seq,'r')
    record = SeqIO.parse(seq_handle,"fasta").next()
    
    insert_node_name = record.id
    
    # inserts the new sequences at the node of the highest scoring HMM
    # needs to change to insert sequence above te highest scoring node
    for node in pt.traverse():
        if not node.is_leaf():
            if node.hmm == highest_node:
                ndist = node.dist/2
                new_node = node.add_child(name=insert_node_name, dist=ndist)
                #R = pt.get_midpoint_outgroup()
                #pt.set_outgroup(R)
                break
    
    # outputs a new tree
    print 'insert_sequence: After'
    pt.write(outfile='final.newick')
    print pt
        
    # aligns the new sequence to the the HMM
    tree_dir = './'+base+'_fat_cat'
    hmm_dir = tree_dir+'/hmm/'
    os.system('hmmalign --allcol  --informat FASTA --amino --outformat A2M -o final.aln '+hmm_dir+'node0.hmm '+seq)  
    # makes a new MSA
    os.system('cat '+aln_name+' >> final.aln')

#def optimize_branch_lengths():

#Uses RaxML to optimize branch lenghths
#raxml (options)
# -f  e              : optimize model+branch lengths for given input tree under GAMMA/GAMMAI only 
# -t  <file.newick>  : Specify a user starting tree file name in Newick format
# -m PROTGAMMAJTT"   : specified AA matrix (JTT) + Optimization of substitution rates + GAMMA model of rate
# -s  <file.phylip>  : Specify the name of the alignment data file in PHYLIP format
# -n  <file.newick>  : Specifies the name of the output file
# -w  <dir>          : Name of the working directory where RAxML will write its output files
 

def score_sequence(seq, base):
    # scores a new sequence against all HMMs on the tree .
    # First use hmmpress to prepare the HMM DB from the concatenated hmm file
    # and then use hmmscan to find the node associated with the highest scoring hmm
    hmm_db = base+'_concat.hmm'
    hmmscan_out = base+'.scanout'
    os.system('hmmpress -f '+hmm_db)
    os.system('hmmscan --tblout '+hmmscan_out+' '+hmm_db+' '+seq)
    scanout = open(hmmscan_out,'r')

    highest_node = ''   
    highest_score = -1.0
    for l in scanout:
        if l[0] != '#':
            tokens = l.split()
            # Note: node name reflects path of corresponding HMM file. Need to extract just the node name 
            node = tokens[0].split('/').pop()
            # Using tokens[5], full sequence score; tokens[8] is the domain score
            score = float(tokens[5])
            if score > highest_score:
                highest_score = score
                highest_node = node
    return highest_node
    

def build_hmm_from_tree(base,tree_name,aln_name):
    # Reads tree and corresponding msa and create an MSA & HMM for each internal node.
    
    # Clean up and create directories for node msa and hmm files
    tree_dir = './'+base+'_fat_cat'
    os.system('rm -R '+tree_dir)
    os.system('mkdir '+tree_dir)
    msa_dir = tree_dir+'/msa/'
    hmm_dir = tree_dir+'/hmm/'
    os.system('mkdir '+msa_dir)
    os.system('mkdir '+hmm_dir)
    
    
    # Annotate internal nodes with name of corresponding HMM.
    pt = PhyloTree(tree_name,alignment=aln_name,alg_format="fasta")
    #print pt
    #R = pt.get_midpoint_outgroup()
    #pt.set_outgroup(R)
    i_node = 0
    for node in pt.traverse():
        if not node.is_leaf():
            node_name = 'node'+str(i_node)
            node.add_features(hmm=node_name)
            i_node += 1
            # make msa for node
            msa = open(msa_dir+node_name+'.aln','w')
            for leaf in node.get_leaves():
                msa.write('>'+leaf.name+'\n')
                msa.write(leaf.sequence+'\n')
            msa.close()
            # build HMM for node
            os.system('python build_hmmer3_hmm_from_alignment.py --name '+hmm_dir+node_name+' '+msa_dir+node_name+'.aln')
    #concatenate HMMs into one file for Hmmscan
    os.system('cat '+hmm_dir+'*.hmm > '+tree_name+'_concat.hmm')
    return pt

def main():
    print 'FAT CATs are better than SKINNY CATs'
    treeName = sys.argv[1]
    alnName = sys.argv[2]
    seqName = sys.argv[3]
    base = treeName.strip('.newick')        
    """
    will work on this later
    
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dirName",
                      help="all files are saved in DIR directory", metavar="DIR")
    parser.add_option("-t" , --tree, dest="treeName",
                      help="tree that will have sequence added to it")
    """
    pt = build_hmm_from_tree(base,treeName,alnName)
    highest_node = score_sequence(seqName, base)
    insert_sequence(pt,highest_node, seqName,base, alnName)
    print 'done'

main()