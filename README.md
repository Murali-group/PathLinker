# PathLinker

This is an implementation the PathLinker software, which infers
signaling pathways in protein interaction networks.

### Authors
  * Christopher L. Poirel (chris.poirel@gmail.com)
  * Anna Ritz (annaritz@vt.edu)
  * Nicholas Sharp (nsharp3@vt.edu)
  * T. M. Murali (tmmurali@cs.vt.edu) [**corresponding author**]

### Relevant Publications
  * Pathways on Demand: Automated Reconstruction of Human Signaling
  Networks. Anna Ritz, Christopher L. Poirel, Allison N. Tegge, Nicholas
  Sharp, Allison Powell, Kelsey Simmons, Shiv D. Kale, and T. M. Murali
  Virginia Tech, Blacksburg, VA. Manuscript under review.
  * Poirel, C. L., Rodrigues, R. R., Chen, K. C., Tyson, J. J., &
  Murali, T. M. (2013). Top-down network analysis to drive bottom-up
  modeling of physiological processes. Journal of Computational Biology,
  20(5), 409-418.

### Overview 

  Consider a protein-protein interaction network as a
  directed, unweighted graph. Given a set of receptors and a set of
  transcription factors in the network which begin and end a signaling
  pathway, we wish to find interactions which are likely to be
  intermediate steps in the pathway. PathLinker ranks interactions in
  the network with respect to such a query.

PathLinker uses two primary steps to produce this ranking:
  * It runs the PageRank algorithm to find visitation probabilities for
  each node. The PageRank uses a modified teleportation weight, so
  teleporting walkers always return to the receptor nodes. These node
  visitation probabilities are used to compute edge traversal
  probabilities, which are then used as weights in the graph for the
  next step.
  * It computes the *k*-shortest simple paths in the network from any
  receptor to any transcription factor. This is done using a new
  modification of Yen's algorithm which allows fast computation for very
  large *k* values. The interactions (edges) in the network are
  ultimately ranked by the index of the first path in which they appear.

See the publications referenced above for a formal description of the
method.

### What's included
  * **PathLinker.py** An end-to-end implementation of the PathLinker
  algorithm. Given a network along with a set or receptors and a set of
  transcription factors, PathLinker outputs a ranking of edges.
  * **PageRank.py** An implementation of the PageRank algorithm. Used as
  a component of PathLinker, but can also be run as a standalone tool.
  Given a network, PageRank iteratively computes the node visitation
  probability by a teleporting random walker.
  * **ksp_Astar.py** An implementation of Yen's *k*-shortest simple
  paths algorithm augmented to use A\* for additional speedup. Used as a
  component of PathLinker, but can also be run as a standalone tool.
  Given a weighted, directed network, the algorithm computes the *k*
  shortest simple paths betwen a source and target node. Here *simple*
  means that there are no repeated nodes in a path.

Run any of the programs with the --help option for full documentation.
The /example directory contains a sample usage of these programs.

### Dependencies Python packages:
  * NetworkX

The code is valid in both Python 2 and Python 3.

Paths of equal length may be returned and processed in a different
order according to your version of Python and NetworkX. PathLinker
defines an equal ranking in this case; the results should be interpreted
as equivalent, although the output is not textually identical.

### License

GNU GPLv3

This is research code, and represents a prototype implementation of a
new idea, not a stable release of extensively verified software. We
encourage experimentation and collaboration with this work, don't
hesitate to contact us! Please cite the relevant publications listed
above if this code contributes to any research work.

