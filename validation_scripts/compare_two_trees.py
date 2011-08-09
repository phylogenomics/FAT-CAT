#!/usr/bin/env python
'''
Purpose - To compare two trees (newick format) and return the distance between them.

Input - -f Paths to the files containing the two trees

         Four optional arguments are available.
         -r gives the topological difference,
         i.e. the sum of the number of splits found in one tree, but not in the other
         -d gives the weighted robinson-foulds distance (this is the default),
         i.e. the sum of the absolute differences in branch lengths for equivalent splits between two trees
         Branch length for a missing split taken to be 0.0.
         -g The geodesic distance between the two trees as calculated by GeoMeTree.
         -i tokens to be ignored
         -t Run unittests on this file.
         
Output - distance or metric
Please run compare_two_trees_test.py for sample test.
If the sample test works, but not this file, then the problem is with the tree formatting.
'''

import os
import sys
import optparse
import subprocess
import dendropy
from ete2 import Tree
    
taxa = dendropy.TaxonSet()

def option_parser():
    parser = optparse.OptionParser(description="Compare two trees and return the distance between them. \
    The trees should be newick-formatted.")
    parser.add_option('-r', help='Calculates the Robinson-foulds metric', dest='r',
                      default=False, action='store_true')
    parser.add_option('-d', help='Calculates the weighted Robinson-foulds distance', dest='d',
                      default=False, action='store_true')
    parser.add_option('-g', help='Calculates the geodesic distance for all splits using GeoMeTree', dest='g',
                      default=False, action='store_true')
    parser.add_option('-f', help='Paths to the two tree files that need to be compared.', dest='treeFiles',
                      action='store', nargs=2, type='string', metavar='<File1> <File2>')
    parser.add_option('-l', help='Compare one leaf node distance between two trees', dest='l',
                      action='store', type='string', metavar='<leafName>')
    parser.add_option('-i',
                      help='Ignore one or more tokens if present in the file string. \
                      Add these as a continuous string e.g. /+, ',
                      dest='ignore_tokens',
                      action='callback', type='string', metavar='e.g. /', callback = list_callback)
    parser.add_option('-t', help='Run unit tests', dest='t',
                      default=False, action='store_true')
    if len(sys.argv)==1:
        parser.print_help()
    return parser.parse_args()

def list_callback(option, opt, value, parser):
    '''
    Getting a list of values from the argument. Code from
    http://stackoverflow.com/questions/392041/python-optparse-list
    '''
    setattr(parser.values, option.dest, list(value))
    
def get_tree_list(files_tuple, ignore_tokens = None):
    '''
    Makes a ree list from the files tuple. Characters to be ignored are passed here.
    '''
    tree_list = [open(x).read() for x in files_tuple]
    tree_list = [create_tree_from_string(x, ignore_tokens) for x in tree_list]
    return tree_list

def create_tree_from_string(treeString, ignore_tokens = None):
    '''
    Creates a dendropy Tree object after replacing the ignored tokens with |
    '''
    if ignore_tokens:
        for token in ignore_tokens:
            treeString = treeString.replace(token, '|')
    return dendropy.Tree.get_from_string(treeString, "newick",
                                         tree_offset=0,
                                         taxon_set=taxa)
    
def robinson_foulds_metric(tree1, tree2):
    '''
    Topological difference by option -r
    This is the base Robinson-foulds metric that is given by dendropy.treecalc.symmetric_difference
    '''
    return dendropy.treecalc.symmetric_difference(tree1, tree2)

def robinson_foulds_distance(tree1, tree2):
    '''
    Distance difference by option -d
    Edge length based weighted Robinson-foulds distance given by dendropy.treecalc.robinson_foulds_distance
    '''
    return dendropy.treecalc.robinson_foulds_distance(tree1, tree2)

def calculate_rf_metric(files_tuple, ignore_tokens = None):
    '''
    Main function called to calculate the robinson-foulds metric
    '''
    tree_list = get_tree_list(files_tuple, ignore_tokens)
    return robinson_foulds_metric(tree_list[0], tree_list[1])

def calculate_rf_distance(files_tuple, ignore_tokens = None):
    '''
    Main function called to calculate the weighted robinson-foulds distance.
    '''
    tree_list = get_tree_list(files_tuple, ignore_tokens)
    return robinson_foulds_distance(tree_list[0], tree_list[1])

def calculate_geodesic_distance(files_tuple, ignore_tokens = None):
    '''
    Distance difference by option -g
    Function using GeoMetree to calculate geodesic distance. Returns a list containing 
    '''
    import os, re
    #Open the two tree files and combine them into one
    tree_1 = open(files_tuple[0]).read()
    tree_2 = open(files_tuple[1]).read()
    combined_file = '%s\n%s' % (tree_1, tree_2)
    f = open('combined.newick', 'w'); f.write(combined_file); f.close()
    #Call GeoMeTree, get the all splits distance and remove intermediate files.
    subprocess.call(['python', 'GeoMeTree/GeoMeTree.py', '-f', 'combined.newick'])
    result_file = open('pair_1_2').read()
    geometree_result = open('pair_1_2').read()
    geometree_result = ''.join(geometree_result.split('\n'))
    geom_all = re.compile('Results\s*for.*?all\s*splits.*?Geodesic\s*distance\s*(.*?)Cone')
    tree_vals = geom_all.findall(geometree_result)
    os.remove('combined.newick'); os.remove('pair_1_2')
    
    if len(tree_vals) > 0:
        return tree_vals[0]
    else:
        return None
    
def calculate_distance_of_one_leaf_between_trees(files_tuple, leaf_label, ignore_tokens=None):
    '''
    Calculating the distance of input leaf from other leaves between two trees.
    e.g.
    Given a true tree with leaves A,B,C,D
    if inserted tree is generated from the true tree by first pruning and then inserting leaf A,
    then a distance function can be
    calculated using the following steps.
    Get distance(A,B) in true tree - distance(A,B) in inserted tree (Absolute value).
    Repeat this for (A,C), (A,D).
    Sum all distances obtained and get the error function
    '''
    if not leaf_label: return "Missing the input leaf argument."
    tree_1 = open(files_tuple[0]).read()
    tree_2 = open(files_tuple[1]).read()
    return calculate_distance_of_one_leaf_string_input(tree_1, tree_2, leaf_label)

def calculate_distance_of_one_leaf_string_input(tree1_str, tree2_str, leaf_label):
    '''
    See description for calculate_distance_of_one_leaf_between_trees.
    '''
    calculated_distance = 0
    tree_1 = Tree(tree1_str)
    tree_2 = Tree(tree2_str)
    for node in tree_1.iter_leaf_names():
        difference_in_leaf_distance = abs(tree_1.get_distance(leaf_label, node)
                                       - tree_2.get_distance(leaf_label, node))
        calculated_distance += difference_in_leaf_distance
    return calculated_distance
    

def run_unittests():
    '''
    Call the unit tests from compare_two_trees_test
    '''
    import compare_two_trees_test
    import unittest
    suite = unittest.TestLoader().loadTestsFromTestCase(compare_two_trees_test.KnownTreeComparisons)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    (opt, args) = option_parser()
    if opt.t:
        run_unittests()
        exit()
    if not opt.treeFiles:
        print 'Two tree files are needed to continue. Use option -f'
        exit()
    if opt.r:
        print calculate_rf_metric(opt.treeFiles, opt.ignore_tokens)
    elif opt.g:
        print calculate_geodesic_distance(opt.treeFiles, opt.ignore_tokens)
    elif opt.d:
        print calculate_rf_distance(opt.treeFiles, opt.ignore_tokens)
    else:
        print calculate_distance_of_one_leaf_between_trees(opt.treeFiles, opt.l)
    # tree1_str = '((A:0.2,(B:0.4, C:1.1):0.1):0.5, (D:0.2,(E:0.4, F:1.1):0.1):0.6);'
    # tree2_str = '((D:0.2,(C:0.4, F:1.1):0.1):0.5, (A:0.2,(B:0.4, E:1.1):0.1):0.6);'
    # leaf_label = 'C'
    # print calculate_distance_of_one_leaf_string_input(tree1_str, tree2_str, leaf_label)
        
