"""
PathLinker finds high-scoring paths that connect a source node set to a
target node set in a large, optionally weighted directed network.  To
reconstruct signaling pathways, PathLinker returns high-scoring paths
from any receptor protein (source node) to any transcriptional regulator
(target node) in a protein-protein interactome.

See the following paper for more details:
Pathways on Demand: Automated Reconstruction of Human Signaling Networks
Anna Ritz, Christopher L. Poirel, Allison N. Tegge, Nicholas Sharp, 
Allison Powell, Kelsey Simmons, Shiv D. Kale, and T. M. Murali
Virginia Tech, Blacksburg, VA
Manuscript under review.

This code is authored by:
Mitch Wagner: mitchw94@vt.edu
Nicholas Sharp: nsharp3@vt.edu
Anna Ritz: annaritz@vt.edu
Christopher L. Poirel: chris.poirel@gmail.com
T. M. Murali: tmmurali@cs.vt.edu

"""

import types
from math import log

import networkx as nx


# Modifies the structure of the graph by removing all edges entering
# sources. These edges will never contribute to a
# path, according to the PathLinker formulation.
def modifyGraphForKSP_removeEdgesToSources(net, sources):

    # We will never use edges leaving a target or entering a source, since
    # we do not allow start/end nodes to be internal to any path.
    for u,v in net.edges():
        if not net.has_edge(u, v):
            continue
        # remove edges coming into sources
        elif v in sources:
            net.remove_edge(u,v)
    return

# Modifies the structure of the graph by removing all edges 
# exiting targets. These edges will never contribute to a
# path, according to the PathLinker formulation.
def modifyGraphForKSP_removeEdgesFromTargets(net, targets):

    # We will never use edges leaving a target or entering a source, since
    # we do not allow start/end nodes to be internal to any path.
    for u,v in net.edges():
        if not net.has_edge(u, v):
            continue
        # remove edges leaving targets
        elif u in targets:
            net.remove_edge(u,v)
    return

# Modifies the structure of the graph by adding a supersource with an
# edge to every source, and a supersink with an edge from every target.
# These artificial edges are given weight 'weightForArtificialEdges',
# which should correspond to a "free" edge in the current interpretation
# of the graph.
def modifyGraphForKSP_addSuperSourceSink(net, sources, targets, weightForArtificialEdges=0):

    # Add a supersource and supersink. Shortest paths from source to sink are the same
    # as shortest paths from "sources" to "targets".
    for s in sources:
        net.add_edge('source', s, weight=1)
        net.edge['source'][s]['ksp_weight'] = weightForArtificialEdges
    for t in targets:
        net.add_edge(t, 'sink', weight=1)
        net.edge[t]['sink']['ksp_weight'] = weightForArtificialEdges
    return

# Applies a user specified edge penalty to each edge.
# This weight penalizes the score of every path by a factor equal
# to (the number of edges in the path)^(this factor).
#
# This was previously done in the logTransformEdgeWeights method
# with a parameter weight=(sum of all edge weights). In the "standard"
# PathLinker case, this was necessary to account for the probability that
# is lost when edges are removed in modifyGraphForKSP_removeEdges(), along
# with probability lost to zero degree nodes in the edge flux calculation.
def applyEdgePenalty(net, weight):

    if weight == 1.0:
        return
     
    for u,v in net.edges():
        w = net.edge[u][v]['ksp_weight']/weight
        net.edge[u][v]['ksp_weight'] = w
    return

    
# Apply a negative logarithmic transformation to edge weights,
# converting multiplicative values (where higher is better) to additive
# costs (where lower is better).
#
# Before the transformation, weights are normalized to sum to one,
# supporting an interpretation as probabilities.
#
# If the weights in the input graph correspond to probabilities,
# shortest paths in the output graph are maximum-probability paths in
# the input graph.
def logTransformEdgeWeights(net):

    for u,v in net.edges():
        w = -log(max([0.000000001, net.edge[u][v]['ksp_weight']]))/log(10)
        net.edge[u][v]['ksp_weight'] = w
    return


# "Undoes" the logarithmic transform to the path lengths, converting
# the path lengths back into terms of the original edge weights
def undoLogTransformPathLengths(paths):

    new_paths_list = []

    # Reconstructs the path list with each edge distance un-log transformed.
    # We build a new list because tuples are unmodifiable.
    for path in paths:
        new_path = [(x[0], 10 ** (-1 * x[1])) for x in path]
        new_paths_list.append(new_path)
    return new_paths_list


# Given a probability distribution over the nodes, calculate the
# probability "flowing" though the outgoing edges of every node. Used
# to assign edge weights after PageRank-ing nodes.
def calculateFluxEdgeWeights(net, nodeWeights):
    
    # the flux score for and edge (u,v) is f_uv = (w_uv p_u)/d_u where
    # w_uv is the weight of the edge (u,v), p_u is the normalized visitation
    # probability (or node score) for u, and d_u is the weighted outdegree of node u.

    # assign EdgeFlux scores to the edges
    for u,v in net.edges():
        w = nodeWeights[u] * net[u][v]['weight']/net.out_degree(u, 'weight')
        net.edge[u][v]['ksp_weight'] = w
    return


# Print the edges in the k shortest paths graph computed by PathLinker.
# This creates a tab-delimited file with one edge per line with three columns:
# tail, head, and KSP index.
# Here, 'ksp index' indicates the index of the first shortest path in which the edge is used.
def printKSPGraph(f, graph):

    outf = open(f, 'w')
    outf.write('#tail\thead\tKSP index\n')
    edges = graph.edges(data=True)

    # Print in increasing order of KSP identifier.
    for e in sorted(edges, key=lambda x: x[2]['ksp_id']):
        t, h, attr_dict = e
        ksp = attr_dict['ksp_id']
        outf.write('%s\t%s\t%d\n' %(t, h, ksp))
    outf.close()
    return


# Print the k shortest paths in order.
# This creates a tab-delimited file with three columns: the number of
# the path, the length of the path (sum of weights), and the sequence of
# nodes in the path. 
def printKSPPaths(f, paths):

    outf = open(f, 'w')
    outf.write('#KSP\tpath_length\tpath\n')

    for k,path in enumerate(paths, 1):
        pathNodes = [n for n,w in path]
        length = path[-1][1]
        outf.write('%d\t%0.5e\t%s\n' %(k, length, '|'.join(pathNodes[1:-1]) ))
    outf.close()
    return

# Print the edges with the flux weight.
# Sort by decreasing of flux weight.
def printEdgeFluxes(f, graph):

    outf = open(f, 'w')
    outf.write('#tail\thead\tedge_flux\n')
    edges = graph.edges(data=True)

    # Print in decreasing  flux weight
    for e in sorted(edges, key=lambda x: x[2]['ksp_weight'], reverse=True):
        t, h, attr_dict = e
        w = attr_dict['ksp_weight']
        outf.write('%s\t%s\t%0.5e\n' %(t, h, w))
    outf.close()
    return

def readNetworkFile(network_file, pagerank=False):
    
    net = nx.DiGraph()

    # Read the network file
    # infile = open(network_file, 'r')
    for line in network_file:
        items = [x.strip() for x in line.rstrip().split('\t')]

        # Skip empty lines or those beginning with '#' comments
        if line=='\n':
            continue
        if line[0]=='#':
            continue

        id1 = items[0]
        id2 = items[1]

        # Ignore self-edges
        if id1==id2:
            continue

        # Possibly use an edge weight
        eWeight = 1

        if (not pagerank):
            if(len(items) > 2):
                eWeight = float(items[2])

            else:
                raise Exception(
                    "ERROR: All edges must have a weight "
                    "unless --PageRank is used. Edge (%s --> %s) does "
                    "not have a weight entry." % (id1, id2))
                
        # Assign the weight. Note in the PageRank case, "weight" is
        # interpreted as running PageRank and edgeflux on a weighted
        # graph. 
        net.add_edge(id1, id2, ksp_weight=eWeight, weight=eWeight)

    return net
