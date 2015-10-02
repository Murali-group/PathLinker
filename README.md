# PathLinker

This software implementats the PathLinker algorithm for reconstructing
signaling pathways from protein interaction networks.

![An image of a signaling pathway](pathway.jpg)

### Authors
  * Christopher L. Poirel (chris.poirel@gmail.com)
  * Anna Ritz (annaritz@vt.edu)
  * Nicholas Sharp (nsharp3@vt.edu)
  * T. M. Murali (tmmurali@cs.vt.edu) [**corresponding author**]

### Relevant Publications
  * Pathways on Demand: Automated Reconstruction of Human Signaling
  Networks. Anna Ritz, Christopher L. Poirel, Allison N. Tegge, Nicholas
  Sharp, Allison Powell, Kelsey Simmons, Shiv D. Kale, and T. M. Murali
  Virginia Tech, Blacksburg, VA. *Manuscript under review*.
  * Poirel, C. L., Rodrigues, R. R., Chen, K. C., Tyson, J. J., &
  Murali, T. M. (2013). Top-down network analysis to drive bottom-up
  modeling of physiological processes. *Journal of Computational Biology*,
  20(5), 409-418.

### Overview 

  PathLinker is a method for reconstructing signaling pathways from
  protein interaction networks. Consider a protein-protein interaction
  network represented as a directed, unweighted graph. The network may
  contain physical, signaling, and regulatory interactions between pairs
  of proteins. Given a query composed of a set of receptors and a set of
  transcription factors in this network that represent the "sources" and
  "targets" of some specific signaling pathway (e.g., the Wnt pathway),
  PathLinker attemps to compute a high quality reconstruction of this
  pathway. Specifically, PathLinker computes a subnetwork that connects
  the receptors in the query to the transcription factors in the query.
  In its output, PathLinker ranks the interactions in the network with
  respect to the query.

The primary algorithmic component of PathLinker is computing the 
*k*-shortest simple paths in the network from any receptor to any
transcription factor. This is accomplished with a novel modification of 
Yen's algorithm which allows fast computation for very large *k* 
values, e.g., 10,000. PathLinker ranks each interaction in the 
network by the index of the first path in which it appears.

PathLinker supports several different options for weighting the network.
If the input graph has weights, these are interepted as multiplicative
probabilities by default, and the most probable paths are computed after
a log-transform. A flag is provided to interpret edge weights as additive
costs instead. If the input graph does not have weights, PathLinker can
run the PageRank algorithm to generate meaningful probablistic weights
from visitation probabilities.

See the publications referenced above for a formal description of the
method.

### What's included
  * **PathLinker.py** An end-to-end implementation of the PathLinker
  algorithm. Given a network along with a set or receptors and a set of
  transcription factors, PathLinker outputs a ranking of edges.
  * **PageRank.py** An implementation of the PageRank algorithm. Used as
  a component of PathLinker, but can also be run as a standalone tool.
  Given a weighted, directed network, PageRank iteratively computes the
  node visitation probability by a teleporting random walker.
  * **ksp_Astar.py** An implementation of Yen's *k*-shortest simple
  paths algorithm augmented to use the A\* algorithm for additional
  speedup. Used as a component of PathLinker, but can also be run as a
  standalone tool.  Given a weighted, directed network, the algorithm
  computes the *k* shortest simple paths betwen a source and target
  node. Here a path is *simple* if there are no repeated nodes in the
  path.

Run any of the programs with the --help option for full documentation.
The /example directory contains a sample usage of these programs.
For all input files, lines beginning with a hash (#)
    are considered to be comments and ignored.

### Dependencies 

Python packages:
  * NetworkX

The code is valid in both Python 2 and Python 3.

Paths of equal length may be returned and processed in a different
order according to your version of Python and NetworkX. The results should be interpreted
as equivalent, although the output is not textually identical.

### License

GNU GPLv3

### How to Cite PathLinker

We will be very glad to hear from you if you use PathLinker in your work. If you publish a paper that uses PathLinker, please cite

Poirel, C. L., Rodrigues, R. R., Chen, K. C., Tyson, J. J., & Murali, T. M. (2013). Top-down network analysis to drive bottom-up modeling of physiological processes. *Journal of Computational Biology*, 20(5), 409-418.

We encourage experimentation with PathLinker. Please don't hesitate to contact us if you would like to collaborate! 
