#!/usr/bin/env python

'''
FAT_CAT (VERSION 1)
Takes a tree (T-), msa (MSA-) and sequence (S) as input.
Builds hmms at every interior node of T-.
Scores S agianst every hmm
Inserts S above the node corresponding to the highest scoring hmm.
Optmizes the branch lengths of the new tree (T- + S)
Aligns S to the highest scoing hmm = S*
Concatenates S* to the root MSA ie MSA- + S*
Outputs final.newick = (T- + S)
        final.aln = (MSA- + S*) 

'''
import os
import sys
import shutil
import optparse
import subprocess
from subprocess import call, check_call
import ete2
from ete2 import PhyloTree
from Bio import SeqIO
from Bio import AlignIO

#Global variables


def option_parser():
    '''
    Options -t for tree file, -a for alignment file, -s for sequence and -o for output directory.
    Called in main()
    '''
    description = 'Fast Approximate Tree classifiCATion. FAT CATs are better than SKINNY CATs.'
    usage = "usage: %prog [options] -t treeFile -m msaFile -s seqFile -o dirName"
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option('-t', help='Path to the tree file (newick-formatted)', dest='t',
                      action='store', metavar="test.newick")
    parser.add_option('-a', help='Path to the alignment msa file (FastA format)', dest='a',
                      action='store', metavar="test.msa")
    parser.add_option('-s', help='Path to the sequence file (unaligned FastA format)', dest='s',
                      action='store', metavar="test.seq")
    parser.add_option('-o', help='Output directory', dest='o',
                      action='store', metavar="testDir")
    if len(sys.argv)==1:
        parser.print_help()
    return parser

def insert_sequence(pt,highest_node,seq, base,aln_name,final_dir):
    '''
    Inserts sequence above the highest scoring node
    '''
    print 'insert_sequence: Before'
    seq_handle = open(seq,'r')
    record = SeqIO.parse(seq_handle,"fasta").next()
    
    insert_node_name = record.id
    
    # insert sequence above te highest scoring node
    for node in pt.traverse():
        if not node.is_leaf():
            if node.hmm == highest_node:
                children = node.get_children()
                for child in children:
                    node.remove_child(child)
                intermediate_node = node.add_child(name='NoName')
                for child in children:
                    intermediate_node.add_child(child)
                new_node = node.add_child(name=insert_node_name)
                break
    
    # outputs a new tree
    print 'insert_sequence: After'
    pt.write(outfile='%sfinal_inserted.newick' % final_dir)
    print pt
        
    # aligns the new sequence to the the HMM
    tree_dir = './%s_fat_cat' % base
    hmm_dir = '%s/hmm/' % tree_dir
    check_call(['hmmalign', '--allcol', '--informat', 'FASTA', '--amino', '--outformat', 'A2M',
          '--trim', '-o', '%sfinal_inserted.aln' % final_dir,
          '%s%s.hmm' % (hmm_dir, highest_node),
          seq])
    os.system('cat %s >> %sfinal_inserted.aln' % (aln_name, final_dir))

def optimize_branch_lengths_fasttree(final_dir):
    '''
    Uses FastTree to optimize branch lengths
    FastTree options:
     -nome -mllen with -intree to optimize branch lengths for a fixed topology
    '''
    input_aln = '%sfinal_inserted.aln' % final_dir
    input_tree = '%sfinal_inserted.newick' % final_dir
    output_tree = "%sfinal_FastTree_optimized.newick" % final_dir
    os.system('FastTree -nome -mllen -intree %s %s > %s' % (input_tree, input_aln, output_tree))
    
def optimize_branch_lengths_raxml(final_dir,base,msa_dir,tree_dir):
    '''
    Uses RaxML to optimize branch lenghths
    raxml (options)
     -f  e              : optimize model+branch lengths for given input tree under GAMMA/GAMMAI only 
     -t  <file.newick>  : Specify a user starting tree file name in Newick format
     -m PROTGAMMAJTT"   : specified AA matrix (JTT) + Optimization of substitution rates + GAMMA model of rate
     -s  <file.phylip>  : Specify the name of the alignment data file in PHYLIP format
     -n  <file.newick>  : Specifies the name of the output file
     -w  <dir>          : Name of the working directory where RAxML will write its output files
    '''
    input_aln = '%sfinal_inserted.aln' % final_dir
    input_tree = '%sfinal_inserted.newick' % final_dir
    
    leaf_map = map_files(input_aln,input_tree,base,msa_dir,tree_dir,final_dir)
    i = open('%sfinal_mapped.aln' %  final_dir)
    o = open('%sfinal_mapped.aln.phylip' % final_dir,'w')
    AlignIO.write(AlignIO.parse(i,'fasta'),o,'phylip')
    i.close()
    o.close()
    
    # need to write to final directory need to use -w option correctly???
    os.system('raxmlHPC -f e -t '+final_dir+'final_mapped.newick'+'  -m PROTGAMMAJTT -s '+final_dir+'final_mapped.aln.phylip -n optimized.newick')
    #os.system('raxmlHPC -g '+final_dir+'final_mapped.newick'+'  -m PROTGAMMAJTT -s '+final_dir+'final_mapped.aln.phylip -n minus_g.newick')
    
def score_sequence(seq, base,hmm_dir):
    '''
    scores a new sequence against all HMMs on the tree .
    First use hmmpress to prepare the HMM DB from the concatenated hmm file
    and then use hmmscan to find the node associated with the highest scoring hmm
    '''
    hmm_db = '%s%s_concat.hmm' % (hmm_dir, base)
    hmmscan_out = '%s%s.scanout' % (hmm_dir, base)
    check_call(['hmmpress', '-f', hmm_db])
    check_call(['hmmscan', '--tblout', hmmscan_out, hmm_db, seq])
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
    

def build_hmm_from_tree(base,tree_name,aln_name,msa_dir,hmm_dir):
    '''
    Reads tree and corresponding msa and create an MSA & HMM for each internal node.
    '''
    
    # Annotate internal nodes with name of corresponding HMM.
    pt = PhyloTree(tree_name,alignment=aln_name,alg_format="fasta")
    i_node = 0
    for node in pt.traverse():
        if not node.is_leaf():
            node_name = 'node%s' % (str(i_node))
            node.add_features(hmm=node_name)
            i_node += 1
            
            # make msa for node
            msa_string = []
            for leaf in node.iter_leaves():
                msa_string.append(">%s" % leaf.name)
                msa_string.append(str(leaf.sequence))
            msa_string = '\n'.join(msa_string)
            msa = open('%s%s.aln' % (msa_dir, node_name),'w'); msa.write(msa_string); msa.close()
            
            # build HMM for node
            check_call(['python', 'build_hmmer3_hmm_from_alignment.py', '--name',
                     '%s%s' % (hmm_dir, node_name),
                     '%s%s.aln' % (msa_dir, node_name)])
                     
    #concatenate HMMs into one file for Hmmscan
    os.system('cat %s*.hmm > %s%s_concat.hmm' %
              (hmm_dir, hmm_dir, base))
    return pt

def map_files(aln_name,tree_name,base,msa_dir,tree_dir,final_dir):
    '''
    Takes the input Tree  msa and maps their names to shorter
    names. The mapping is recorded in a dictionary. This is required for
    branch optimization
    
    '''
    leaf_map = {}
    map_file = open(msa_dir+'msa.map','w')
    
    # Create dictionary of mappings from short name to long name
    mapped_msa_name = final_dir+'final_mapped.aln'
    mapped_msa = open(mapped_msa_name,'w')
    leaf_id = 0
    # Create mapped MSA
    for record in SeqIO.parse(open(aln_name,"rU"),'fasta'):   
        short_name = 'seq_'+str(leaf_id)
        leaf_map[record.id] = short_name
        leaf_id += 1
        map_file.write(record.id+' '+short_name+'\n')
        mapped_msa.write('>'+short_name+'\n'+str(record.seq)+'\n')
        
    # Create mapped Newick
    tree = open(tree_name,'r')
    newick_string = tree.read()
    for leaf_name in leaf_map:
        newick_string = newick_string.replace(leaf_name,leaf_map[leaf_name])
    mapped_tree_name = final_dir+'final_mapped.newick'
    mapped_tree = open(mapped_tree_name,'w')
    mapped_tree.write(newick_string)
    tree.close()
    mapped_tree.close()
    return leaf_map
    
def unmap_files(msa_dir,final_dir):
    msa_map_file = open(msa_dir+'msa.map','r')
    msa_map = []
    for line in msa_map_file:
        line = line.strip()
        map_fields = line.split(' ')
        map_original = map_fields[0]
        map_new = map_fields[1]
        msa_map.append((map_new,map_original))
    
    # unmap Newick
    mapped_newick_file = open('RAxML_result.optimized.newick','r')
    #mapped_newick_file = open('RAxML_result.minus_g.newick','r')
    newick_string = mapped_newick_file.read()
    
    for map in reversed(msa_map):
        newick_string = newick_string.replace(map[0],map[1])
        
    mapped_newick_file.close()
    
    unmapped_newick_file = open(final_dir+'final.newick','w')
    #unmapped_newick_file = open(final_dir+'final_g.newick','w')
    unmapped_newick_file.write(newick_string)
    unmapped_newick_file.close()
    
    # unmap aln
    mapped_aln_file = open(final_dir+'final_mapped.aln','r')
    aln_string = mapped_aln_file.read()
    
    for map in reversed(msa_map):
        aln_string = aln_string.replace(map[0],map[1])
        
    mapped_aln_file.close()
    
    unmapped_aln_file = open(final_dir+'final.aln','w')
    #unmapped_aln_file = open(final_dir+'final_g.aln','w')
    unmapped_aln_file.write(aln_string)
    unmapped_aln_file.close()
    
        
    
def main():
    # Get input arguments
    parser = option_parser()
    (opt, args) = parser.parse_args()
    tree_name = opt.t
    aln_name = opt.a
    seq_name = opt.s
    base = opt.o
    if not (tree_name and aln_name and seq_name):
        parser.print_help()
        print 'Error: Not all inputs have been provided'
        return
    # Clean up and create directories for node msa and hmm files
    main_dir = './%s_fat_cat' % base
    try:
        shutil.rmtree(main_dir)
    except:
        pass
    dir_structure = ['%s/%s/' % (main_dir,d) for d in 'msa', 'hmm', 'tree', 'final']
    msa_dir, hmm_dir, tree_dir, final_dir = dir_structure
    [os.makedirs(d) for d in dir_structure]
    
    
    pt = build_hmm_from_tree(base,tree_name,aln_name,msa_dir,hmm_dir)
    highest_node = score_sequence(seq_name, base,hmm_dir) 
    insert_sequence(pt,highest_node, seq_name,base, aln_name,final_dir)
    optimize_branch_lengths_fasttree(final_dir)
    #optimize_branch_lengths_raxml(final_dir,base,msa_dir,tree_dir)
    #unmap_files(msa_dir,final_dir)
    print highest_node
    print 'done'

if __name__ == "__main__":
    main()
