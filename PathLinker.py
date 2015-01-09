"""
PathLinker finds high-scoring paths that connect a source node set to a
target node set in a large, optionally weighted directed network.  To
reconstruct signaling pathways, PathLinker returns high-scoring paths
from any receptor protein (source node) to any transcriptional regulator
(target node) in a protein-protein interactome.

See the following paper for more details:
Pathways on Demand: Automated Reconstruction of Human Signaling Networks
Anna Ritz, Christopher L. Poirel, Allison N. Tegge, Nicholas Sharp, Allison Powell, Kelsey Simmons, Shiv D. Kale, and T. M. Murali
Virginia Tech, Blacksburg, VA
Manuscript under review.

This code is authored by:
Nicholas Sharp: nsharp3@vt.edu
Anna Ritz: annaritz@vt.edu
Christopher L. Poirel: chris.poirel@gmail.com
T. M. Murali: tmmurali@cs.vt.edu

"""

import sys
from optparse import OptionParser, OptionGroup
import types
from math import log

import networkx as nx

# local modules
import ksp_Astar as ksp
from PageRank import pagerank

# Prepare the network for the KSP algorithm by: 
#   - Modifying the weights along the edges according to the pagerank probability flux
#   - Removing outgoing edges from the TRs and incoming edges to the receptors. 
#   - Creating a 'supersource' node with edges to each of the receptors
#   - Creating a 'supersink' node with edges from each of the TRs
#
# The weights to be used for KSP will be stored in the 'ksp_weight' attribute of the 
# Networkx graph.
def prepareNetForKSP(net, sources, targets, nodeWeights):

    # the flux score for and edge (u,v) is f_uv = (w_uv p_u)/d_u where
    # w_uv is the weight of the edge (u,v), p_u is the normalized visitation
    # probability (or node score) for u, and d_u is the weighted outdegree of node u.
    currWeight = {}

    # we will never use edges leaving a target or entering a source, since
    # we do not allow start/end nodes to be internal to any path.
    for u,v in net.edges():
        if not net.has_edge(u, v):
            continue
        # remove edges leaving targets
        elif u in targets:
            net.remove_edge(u,v)
        # remove edges coming into sources
        elif v in sources:
            net.remove_edge(u,v)

    # assign EdgeFlux scores to the edges
    sumWeight = 0
    for u,v in net.edges():
        currWeight[(u,v)] = nodeWeights[u] * net[u][v]['weight']/net.out_degree(u, 'weight')
        sumWeight += currWeight[(u,v)] # Normalize to account for probability lost from removed edges

    for u,v in net.edges():
        w = -log(max([0.000000001, currWeight[(u,v)] / sumWeight]))/log(10)
        net.edge[u][v]['ksp_weight'] = w

    # add supersource and supersink. shortest paths from source to sink are the same
    # as shortest paths from "sources" to "targets".
    for s in sources:
        net.add_edge('source', s, weight=1)
        net.edge['source'][s]['ksp_weight'] = 0
    for t in targets:
        net.add_edge(t, 'sink', weight=1)
        net.edge[t]['sink']['ksp_weight'] = 0
    return



# Print the edges in the k shortest paths graph computed by Linker.
# This creates a tab-delimited file with one edge per line with three columns:
# tail, head, ksp
# Here, ksp indicates the first shortest path in which the edge is used.
def printKSPGraph(f, graph):

    outf = open(f, 'w')
    outf.write('#tail\thead\tksp\n')
    edges = graph.edges(data=True)

    # Print in increasing order of ksp identifier.
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
    outf.write('#ksp\tpath_length\tpath\n')

    for k,path in enumerate(paths, 1):
        pathNodes = [n for n,w in path]
        length = path[-1][1]
        outf.write('%d\t%0.5e\t%s\n' %(k, length, '|'.join(pathNodes[1:-1]) ))
    outf.close()
    return


def main(args):
    usage = '''
PathLinker.py [options] NETWORK NODE_TYPES
REQUIRED arguments:
    NETWORK - A tab-delimited file with one directed interaction per line. Each
        line should have at least 2 columns: tail, head. Edges are directed from
        tail->head. This file can optionally have a third column specifying the
        edge weight

    NODE_TYPES - A tab-delimited file denoting nodes as receptors or TRs. The first
        column is the node name, the second is the node type, either 'source'
        (or 'receptor') or 'target' (or 'tr' or 'tf'). Nodes which are neither receptors nor TRs may
        be omitted from this file or may be given a type which is neither 'source'
        nor 'target'.

NOTE: This operates on only the largest (weakly) connected component of the input graph.
'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('-o', '--output', type='string', default='out_', metavar='STR',\
        help='A string to prepend to all output files. (default="out")')

    parser.add_option('', '--write-paths', action='store_true', default=False,\
        help='If given, also output a list of paths found by ksp in addition to the ranked edges.')

    # Random Walk Group
    group = OptionGroup(parser, 'Random Walk Options')

    group.add_option('-q', '--q-param', action='store', type='float', default=0.5,\
        help='The value of q indicates the probability that the random walker teleports back to a source node during the random walk process. (default=0.5)')

    group.add_option('-e', '--epsilon', action='store', type='float', default=0.01,\
        help='A small value used to test for convergence of the iterative implementation of PageRank. (default=0.01)')

    group.add_option('', '--max-iters', action='store', type='int', default=500,\
        help='Maximum number of iterations to run the PageRank algorithm. (default=500)')

    parser.add_option_group(group)

    # k shortest paths Group
    group = OptionGroup(parser, 'k Shortest Paths Options')

    group.add_option('-k', '--k-param', type='int', default=100,\
        help='The number of shortest paths to find. (default=100)')

    parser.add_option_group(group)


    # parse the command line arguments
    (opts, args) = parser.parse_args()

    # get the required arguments
    num_req_args = 2
    if len(args)!=num_req_args:
        parser.print_help()
        sys.exit('\nERROR: PathLinker.py requires %d positional arguments, %d given.' %(num_req_args, len(args)))
    NETWORK_FILE = args[0]
    NODE_VALUES_FILE = args[1]

    ## Read the network from file
    net = nx.DiGraph()

    # Read the network file
    print('\nReading the network from %s' %(NETWORK_FILE))
    infile = open(NETWORK_FILE, 'r')
    for line in infile:
        items = [x.strip() for x in line.rstrip().split('\t')]

        # Skip empty lines or those beginning with '#' comments
        if line=='':
            continue
        if line[0]=='#':
            continue

        id1 = items[0]
        id2 = items[1]

        # Ignore self-edges
        if id1==id2:
            continue

        # If no weight is given for the edge, assign it a weight of 1.
        eWeight = 1
        if(len(items) > 2):
            eWeight = float(items[2])

        net.add_edge(id1, id2, weight=eWeight)


    # Operate on only the largest connected component
    print(nx.info(net))
    conn_comps = nx.weakly_connected_component_subgraphs(net)
    
    # This is the only portion of the program which prevents
    # compatability between Python 2 & 3. In 2, this object is a
    # generator, but in 3 it is a list.
    if(isinstance(conn_comps, types.GeneratorType)):
        net = next(conn_comps)
    elif(isinstance(conn_comps, list)):
        net = conn_comps[0]
    else:
        print("Compatability error. Connected components object from NetworkX does not have acceptable type")
        exit(-1)
    

    print("\nLargest weakly connected component:\n"+nx.info(net))

    # Read the sources and targets on which to run pagerank and KSP
    sources = set()
    targets = set()

    # Read the receptors and TRs file
    print("Reading sources and targets from " + NODE_VALUES_FILE)
    for line in open(NODE_VALUES_FILE, 'r').readlines():
        items = [x.strip() for x in line.rstrip().split('\t')]

        # Skip empty lines and lines beginning with '#' comments
        if line=='':
            continue
        if line[0]=='#':
            continue

        if items[1] in ['source', 'receptor']:
            sources.add(items[0])
        elif items[1] in ['target', 'tr', 'tf']:
            targets.add(items[0])

    print('\nRead %d sources and %d targets' % (len(sources), len(targets)))

    # Remove sources and targets that don't appear in the network, and do some sanity checks on sets
    sources = set([s for s in sources if s in net])
    targets = set([t for t in targets if t in net])
    print('\tAfter removing sources and targets that are not in the network: %d sources and %d targets.' %(len(sources), len(targets)))
    if len(sources)==0:
        sys.exit('ERROR: No sources are in the network.')
    if len(targets)==0:
        sys.exit('ERROR: No targets are in the network.')
    if len(sources.intersection(targets))>0:
        sys.exit('ERROR: %d proteins are listed as both a source and target.' %(len(sources.intersection(targets))))


    ## Run PageRank on the network

    PR_PARAMS = {'q' : opts.q_param,\
                 'eps' : opts.epsilon,\
                 'maxIters' : opts.max_iters}

    print('\nRunning PageRank on net.(q=%f)' %(opts.q_param))

    # The initial weights are entirely on the source nodes, so this 
    # corresponds to a random walk that teleports back to the sources.
    weights = dict.fromkeys(sources, 1.0)
    prFinal = pagerank(net, weights, **PR_PARAMS)

    ## Run KSP on the network

    # Weight the edges in the network for KSP computation, as well as
    # removing outgoing edges from targets and incoming edges from
    # sources. Create a supersource and supersink.
    prepareNetForKSP(net, sources, targets, prFinal)
    
    # Run the pathfinding algorithm
    print('\nComputing the k=%d shortest simple paths.' %(opts.k_param))
    paths = ksp.k_shortest_paths_yen(net, 'source', 'sink', opts.k_param, weight='ksp_weight')

    if len(paths)==0:
        sys.exit('\tERROR: Targets are not reachable from the sources.')

    ## Use the results of KSP to rank edges

    # Prepare the k shortest paths for output to flat files
    pathgraph = nx.DiGraph()
    for k,path in enumerate(paths, 1):

        # Process the edges in this path
        edges = []
        for i in range(len(path)-1):
            t = path[i][0]
            h = path[i+1][0]

            # Skip edges that have already been seen in an earlier path
            if pathgraph.has_edge(t, h):
                continue

            # Skip edges that include our artificial supersource or
            # supersink
            if t=='source' or h=='sink':
                continue

            # This is a new edge. Add it to the list and note which k it
            # appeared in.
            else:
                edges.append( (t,h,{'ksp_id':k, 'ksp_weight':net.edge[t][h]['ksp_weight']}) )

        # Add all new, good edges from this path to the network
        pathgraph.add_edges_from(edges)

        # Each node is ranked by the first time it appears in a path.
        # Identify these by check for any nodes in the graph which do
        # not have 'ksp_id' attribute, meaining they were just added
        # from this path.
        for n in pathgraph.nodes():
            if 'ksp_id' not in pathgraph.node[n]:
                pathgraph.node[n]['ksp_id'] = k


    ## Write out the results to file

    # Write a list of all edges encountered, ranked by the path they
    # first appeared in.
    kspGraphOutfile = '%sk_%d-ranked-edges.txt' %(opts.output, opts.k_param)
    printKSPGraph(kspGraphOutfile, pathgraph)
    print('\nKSP results are in "%s"' %(kspGraphOutfile))

    # Write a list of all paths found by the ksp algorithm, if
    # requested.
    if(opts.write_paths):
        kspOutfile = '%sk_%d-paths.txt' %(opts.output, opts.k_param)
        printKSPPaths(kspOutfile, paths)
        print('KSP paths are in "%s"' %(kspOutfile))

    print('\nFinished!')

if __name__=='__main__':
    main(sys.argv)
