'''
Test suite for the compare_two_trees.py file.
Tests situations where two trees are the same, one node has been switched, and two nodes have been switched.

'''

import compare_two_trees
import unittest

class KnownTreeComparisons(unittest.TestCase):
    '''
    Base class for unit tests. Defines several trees for test.
    '''
    orig_tree = '((A:0.2,(B:0.4, C:1.1):0.1):0.5, (D:0.2,(E:0.4, F:1.1):0.1):0.6);'
    same_tree = '((A:0.2,(B:0.4, C:1.1):0.1):0.5, (D:0.2,(E:0.4, F:1.1):0.1):0.6);'
    # one_node_switched_tree should give an rf metric of 4, ie ABC|DEF, BCAD|EF, EBC|DAF, BCED|AF
    # weighted rf distance should be the summation of the weights of these splits i.e., ABC|DEF(1.1),
    # BCAD|EF(0.1), EBC|DAF(1.1), BCED|AF(0.1) and the difference in edge-lengths for A|BCDEF(0.2)
    # and E|ABCDF(0.2) = 2.8
    one_node_switched_tree = '((E:0.2,(B:0.4, C:1.1):0.1):0.5, (D:0.2,(A:0.4, F:1.1):0.1):0.6);'
    # two_nodes_switched_tree should give an rf metric of 6 and a weighted distance of 4.0 (as above)
    two_nodes_switched_tree = '((D:0.2,(C:0.4, F:1.1):0.1):0.5, (A:0.2,(B:0.4, E:1.1):0.1):0.6);'


    # Convert the strings to trees. No tests needed for this one if dendropy is installed correctly.
    orig_tree = compare_two_trees.create_tree_from_string(orig_tree)
    same_tree = compare_two_trees.create_tree_from_string(same_tree)
    one_node_switched_tree = compare_two_trees.create_tree_from_string(one_node_switched_tree)
    two_nodes_switched_tree = compare_two_trees.create_tree_from_string(two_nodes_switched_tree)
    
    def test_original_tree_with_duplicate_tree_rf_metric(self):                          
        """Rf metric of original tree and duplicate tree"""
        result = compare_two_trees.robinson_foulds_metric(self.orig_tree, self.same_tree) 
        self.assertEqual(0, result)

    def test_original_tree_with_duplicate_tree_rf_distance(self):                          
        """Rf weighted distance of original tree and duplicate tree"""
        result = compare_two_trees.robinson_foulds_distance(self.orig_tree, self.same_tree) 
        self.assertEqual(0.0, result)  

    def test_original_tree_with_one_node_switched_tree_rf_metric(self):                          
        """Rf metric of original tree and a duplicate tree with one node switched"""
        result = compare_two_trees.robinson_foulds_metric(self.orig_tree,
                                                                    self.one_node_switched_tree) 
        self.assertEqual(4, result)

    def test_original_tree_with_two_nodes_switched_tree_rf_metric(self):
        """Rf metric of original tree and a duplicate tree with two nodes switched"""
        result = compare_two_trees.robinson_foulds_metric(self.orig_tree,
                                                                    self.two_nodes_switched_tree) 
        self.assertEqual(6, result)
        
    def test_original_tree_with_one_node_switched_tree_rf_distance(self):                          
        """Rf weighted distance of original tree and a duplicate tree with one node switched"""
        result = compare_two_trees.robinson_foulds_distance(self.orig_tree,
                                                                      self.one_node_switched_tree) 
        self.assertEqual("%.2f" % 2.8, "%.2f" % result)

    def test_original_tree_with_two_nodes_switched_tree_rf_distance(self):                          
        """Rf weighted distance of original tree and a duplicate tree with two nodes switched"""
        result = compare_two_trees.robinson_foulds_distance(self.orig_tree,
                                                                      self.two_nodes_switched_tree) 
        self.assertEqual("%.2f" % 4.0, "%.2f" % result)
        

if __name__ == "__main__":
    unittest.main()
