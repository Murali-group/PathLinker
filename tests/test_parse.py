import unittest
from io import StringIO

import parse


class TestParse(unittest.TestCase):


    def setUp(self):
        self.edges_contents = StringIO(
            "#tail\thead\tKSP Index\n"
            "A\tB\t1\n"
            "A\tE\t2\n"
            "E\tF\t2")

        self.paths_contents = StringIO(
            "#KSP\tpath_length\tpath\n"
            "1\t1\tA|B\n"\
            "2\t2\tA|E|F")
   

    def test_parse_ranked_edges(self):
        ranked_edges = parse.parse_ranked_edges(self.edges_contents)

        expected = [set([("A", "B")]), set([("A", "E"), ("E", "F")])]

        self.assertEqual(len(expected), len(ranked_edges))

        for pair in zip(expected, ranked_edges):
            self.assertEqual(pair[0], pair[1])


    def test_parse_ranked_paths(self):
        ranked_paths = parse.parse_ranked_paths(self.paths_contents)

        expected = [("A|B", "1"), ("A|E|F", "2")]

        self.assertEqual(len(ranked_paths), len(expected))

        for pair in zip(expected, ranked_paths):
            self.assertEqual(pair[0], pair[1])


    def test_get_ranked_nodes(self):
        ranked_edges = [set([("A", "B")]), set([("A", "E"), ("E", "F")])]
        ranked_nodes = parse.get_ranked_nodes(ranked_edges)

        expected = [set(["A", "B"]), set(["E", "F"])]

        self.assertEqual(len(ranked_nodes), len(expected))
        
        for pair in zip(expected, ranked_nodes):
            self.assertEqual(pair[0], pair[1])
