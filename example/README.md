# A simple demo of this codebase

This directory contains a small example graph, used to demonstrate the
functionality of this code. There are two primary input files:

  * **sample-in-net.txt** This file contains a description the directed
  graph as a list of weighted edges.
  * **sample-in-nodetypes.txt** This file contains node type information
  to indicate which nodes are sources and which are targets in the
  network.

In the sample network, the source nodes are named "S#" and the target
nodes are named "T#" for clarity. However, these names are completely
arbitrary. Any node names can be used, the program knows that these
nodes are sources/targets because they are indicated as such in the
nodetypes.txt file.

### Running PathLinker

The sample data describes a network with two receptors and three
transcription factors, in addition to 6 other proteins. There are a
total of 20 given interactions in this network. All interactions have an
arbitrary probability weight of 0.5. PathLinker will provide a ranking
for these interactions by how likely that are to appear on a signaling 
pathway between the receptors and transcription factors.

To run PathLinker on this network, use the command:
(assuming the example directory is the working directory)

    python ../PathLinker.py sample-in-net.txt sample-in-nodetypes.txt

This will create a file, `out_k_100-ranked-edges.txt`, which ranks the
interactions (edges) as described above. Some edges are not listed in
this file; PathLinker implicitly ranks these unlisted interactions as
equal, below the explicitly ranked edges.

Adding the option `--write-paths` will generate a second output file
`out_k_100-paths.txt`, which lists the paths used by PathLinker to
determine the ranking. Additional options control the parameters for
PathLinker, as described in the reference materials.

Two sample output files are provided for comparison:
`sample-out-ranked-edges.txt` and `sample-out-paths.txt`. Note that
these files may not exactly match your outputs, because paths of
equal length may be returned and processed in a different order
according to your version of Python and NetworkX. However, recall that
PathLinker defines an equal ranking in this case; the results
should be interpreted as equivalent, although the output is not
textually identical.

### Running PathLinker with PageRank

Suppose you wanted to run PathLinker on that same network, but did not have
edge weights. PathLinker can first run PageRank to generate edge weights.
The file `sample-in-net-noweights.txt` is identical to the first network,
but lacks weights.

To run PathLinker with PageRank on this network, use the command:
(assuming the example directory is the working directory)

    python ../PathLinker.py --PageRank sample-in-net-noweights.txt sample-in-nodetypes.txt

### Running PageRank

To only compute a PageRank on the network, use the command:

    python ../PageRank.py sample-in-net.txt

This will output a file, `pageranks.txt`, which contains a visitation
probability for every node on the network.

To use a teleportation weight file, use the command:

    python ../PageRank.py sample-in-net.txt --tele-weights sample-in-teleprobs.txt

This will output a new `pageranks.txt` for the modified system. Note
that the weights in this input file do not necessarily need to be
normalized, they will be normalized internally. Any nodes omitted from
this file are assumed to have weight 0.

Additional arguments to the program can adjust the algorithm parameters,
as described in the reference materials.

### Running k-shortest-paths

To compute the k shortest simple paths in the network, use the command:

    python ../ksp_Astar.py sample-in-net.txt G E

This will compute the (default k = 200) shortest paths from node G to
node E in the graph. The result is printed in the file `paths.txt`. The
graph only actually contains 4 distinct simple paths from node G to node
E, so only 4 paths are returned.

Additional arguments to the program can adjust the output location and k
value.
