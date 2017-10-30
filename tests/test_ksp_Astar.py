import unittest
import networkx as nx
import ast
from io import StringIO

import ksp_Astar as ksp

# TODO: PathLinker will return results in different orders
# sometimes. For the time being, re-running the tests
# until the pass should be verification enough.

class Test_k_shortest_paths_yen(unittest.TestCase):

    def test_empty_graph(self):
        G = nx.DiGraph()
        paths = ksp.k_shortest_paths_yen(G, '1', '4', 1)
        self.assertEqual([], paths)

    # Test behavior when obtaining a single path
    def test_one_path_k1(self):
        G = nx.DiGraph()
        G.add_edge('1', '2', weight=1)
        G.add_edge('2', '3', weight=1)
        G.add_edge('3', '4', weight=1)
        paths = ksp.k_shortest_paths_yen(G, '1', '4', 1)
        expected = ast.literal_eval(
            "[[('1',0), ('2', 1), ('3', 2), ('4', 3)]]")
        self.assertEqual(paths, expected)

    # Test behavior when obtaining multiple paths where only one exists
    def test_one_path_k10(self):
        G = nx.DiGraph()
        G.add_edge('1', '2', weight=1)
        G.add_edge('2', '3', weight=1)
        G.add_edge('3', '4', weight=1)
        paths = ksp.k_shortest_paths_yen(G, '1', '4', 10)

    def test_disjoint_equal_paths(self):
        G = nx.DiGraph()
        G.add_edge('1', '2', weight=1)
        G.add_edge('2', '3', weight=1)
        G.add_edge('3', '4', weight=1)

        G.add_edge('1', '5', weight=1)
        G.add_edge('5', '6', weight=1)
        G.add_edge('6', '4', weight=1)

        G.add_edge('1', '7', weight=1)
        G.add_edge('7', '8', weight=1)
        G.add_edge('8', '4', weight=1)

        paths = ksp.k_shortest_paths_yen(G, '1', '4', 3)
        expected = ast.literal_eval(
            "[[('1',0), ('2', 1), ('3', 2), ('4', 3)],"
            "[('1',0), ('5', 1), ('6', 2), ('4', 3)],"
            "[('1',0), ('7', 1), ('8', 2), ('4', 3)]]")

        self.assertEqual(paths, expected)

    def test_disjoint_equal_paths_clip(self):
        G = nx.DiGraph()
        G.add_edge('1', '2', weight=1)
        G.add_edge('2', '3', weight=1)
        G.add_edge('3', '4', weight=1)

        G.add_edge('1', '5', weight=1)
        G.add_edge('5', '6', weight=1)
        G.add_edge('6', '4', weight=1)

        G.add_edge('1', '7', weight=1)
        G.add_edge('7', '8', weight=1)
        G.add_edge('8', '4', weight=1)

        paths = ksp.k_shortest_paths_yen(G, '1', '4', 1, clip=True)
        expected = ast.literal_eval(
            "[[('1',0), ('2', 1), ('3', 2), ('4', 3)]]")

        self.assertEqual(paths, expected)

    def test_disjoint_equal_paths_no_clip(self):
        G = nx.DiGraph()
        G.add_edge('1', '2', weight=1)
        G.add_edge('2', '3', weight=1)
        G.add_edge('3', '4', weight=1)

        G.add_edge('1', '5', weight=1)
        G.add_edge('5', '6', weight=1)
        G.add_edge('6', '4', weight=1)

        G.add_edge('1', '7', weight=1)
        G.add_edge('7', '8', weight=1)
        G.add_edge('8', '4', weight=1)

        paths = ksp.k_shortest_paths_yen(G, '1', '4', 1, clip=False)
        expected = ast.literal_eval(
            "[[('1',0), ('2', 1), ('3', 2), ('4', 3)],"
            "[('1',0), ('5', 1), ('6', 2), ('4', 3)],"
            "[('1',0), ('7', 1), ('8', 2), ('4', 3)]]")

        self.assertEqual(paths, expected)

    def test_simple_graph_overlapping_paths(self):
        G = nx.DiGraph()
        G.add_edge('1', '2', weight=1)
        G.add_edge('2', '3', weight=1)
        G.add_edge('3', '4', weight=1)
        G.add_edge('4', '5', weight=1)

        G.add_edge('1', '6', weight=2)
        G.add_edge('6', '4', weight=1.5)
        
        G.add_edge('1', '3', weight=5)
        G.add_edge('1', '4', weight=7)

        paths = ksp.k_shortest_paths_yen(G, '1', '5', 5, clip=True)
        expected = ast.literal_eval(
            "[[('1',0), ('2', 1), ('3', 2), ('4', 3), ('5', 4)],"
            "[('1',0), ('6', 2), ('4', 3.5), ('5', 4.5)],"
            "[('1',0), ('3', 5), ('4', 6), ('5', 7)],"
            "[('1',0), ('4', 7), ('5', 8)]]")
        self.assertEqual(paths, expected)


if __name__ == '__main__':
    unittest.main()
