import unittest
import networkx as nx
from io import StringIO

import PathLinker as pl

class TestReadNetworkFile(unittest.TestCase):

    def test_empty_file(self):
        contents = StringIO(u'')
        net = pl.readNetworkFile(contents, True) 
        self.assertEqual(type(net), nx.DiGraph)
        self.assertEqual(len(net.nodes()), 0)
        self.assertEqual(len(net.edges()), 0)

    def test_small_network_with_weights(self):
        contents = StringIO(u'1\t2\t.5\n1\t3\t.5')
        net = pl.readNetworkFile(contents, True) 
        self.assertEqual(type(net), nx.DiGraph)
        self.assertEqual(len(net.nodes()), 3)
        self.assertEqual(len(net.edges()), 2)

    def test_small_network_empty_lines_are_ignored(self):
        contents = StringIO(u'\n\n\n\n\n1\t2\t.5\n\n1\t3\t.5\n\n\n')
        net = pl.readNetworkFile(contents, True) 
        self.assertEqual(type(net), nx.DiGraph)
        self.assertEqual(len(net.nodes()), 3)
        self.assertEqual(len(net.edges()), 2)

    def test_small_network_commented_lines_are_ignored(self):
        contents = StringIO(u'#comment\n1\t2\t.5\n#comment\n1\t3\t.5')
        net = pl.readNetworkFile(contents, True) 
        self.assertEqual(type(net), nx.DiGraph)
        self.assertEqual(len(net.nodes()), 3)
        self.assertEqual(len(net.edges()), 2)

    def test_self_edges_are_ignored(self):
        contents = StringIO(u'1\t1\t.5\n1\t3\t.5')
        net = pl.readNetworkFile(contents, True) 
        self.assertEqual(type(net), nx.DiGraph)
        self.assertEqual(len(net.nodes()), 2)
        self.assertEqual(len(net.edges()), 1)


    def test_no_weights_pagerank_false_exit(self):
        contents = StringIO(u'1\t1\n1\t3')
        with self.assertRaises(Exception):
            net = pl.readNetworkFile(contents, False) 
        
    def test_no_weights_pagerank_true(self):
        contents = StringIO(u'1\t1\n1\t3')
        net = pl.readNetworkFile(contents, True) 
        self.assertEqual(type(net), nx.DiGraph)
        self.assertEqual(len(net.nodes()), 2)
        self.assertEqual(len(net.edges()), 1)

if __name__ == '__main__':
    unittest.main()
