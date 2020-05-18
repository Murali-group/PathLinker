import sys
import types
import networkx as nx
from optparse import OptionParser, OptionGroup

# Local package imports
import PathLinker as pl 
import ksp_Astar as ksp
from PageRank import pagerank, writePageRankWeights

def main(args):
    usage = '''
run.py [options] NETWORK NODE_TYPES
REQUIRED arguments:
    NETWORK - A tab-delimited file with one directed interaction per
        line. Each line should have at least 2 columns: tail, head. Edges
        are directed from tail to head. This file can have a third column
        specifying the edge weight, which is required unless the --PageRank
        option is used (see --PageRank help for a note on these weights).
        To run PathLinker on an unweighted graph, set all edge weights
        to 1 in the input network.

    NODE_TYPES - A tab-delimited file denoting nodes as receptors or TRs. The first
        column is the node name, the second is the node type, either 'source'
        (or 'receptor') or 'target' (or 'tr' or 'tf'). Nodes which are neither receptors nor TRs may
        be omitted from this file or may be given a type which is neither 'source'
        nor 'target'.
    
'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('-o', '--output', type='string', default='out', metavar='STR',\
        help='A string to prepend to all output files. (default="out")')

    parser.add_option('', '--write-paths', action='store_true', default=False,\
        help='If given, also output a list of paths found by KSP in addition to the ranked edges.')
    
    parser.add_option('', '--no-log-transform', action='store_true', default=False,\
        help='Normally input edge weights are log-transformed. This option disables that step.')

    parser.add_option('', '--largest-connected-component', action='store_true', default=False,\
        help='Run PathLinker on only the largest weakly connected component of the graph. May provide performance speedup.')
        
    parser.add_option('', '--edge-penalty', type='float', default=1.0,\
        help='Factor by which to divide every edge weight. The effect of this option is to penalize the score of every path by a factor equal to (the number of edges in the path)^(weight). (default=1.0)')

    # Random Walk Group
    group = OptionGroup(parser, 'Random Walk Options')

    group.add_option('', '--PageRank', action='store_true', default=False,\
        help='Run the PageRank algorithm to generate edge visitation flux values, which are then used as weights for KSP. A weight column in the network file is not needed if this option is given, as the PageRank visitation fluxes are used for edge weights in KSP. If a weight column is given, these weights are interpreted as a weighted PageRank graph.')

    group.add_option('-q', '--q-param', action='store', type='float', default=0.5,\
        help='The value of q indicates the probability that the random walker teleports back to a source node during the random walk process. (default=0.5)')

    group.add_option('-e', '--epsilon', action='store', type='float', default=0.0001,\
        help='A small value used to test for convergence of the iterative implementation of PageRank. (default=0.0001)')

    group.add_option('', '--max-iters', action='store', type='int', default=500,\
        help='Maximum number of iterations to run the PageRank algorithm. (default=500)')

    parser.add_option_group(group)

    # k shortest paths Group
    group = OptionGroup(parser, 'k Shortest Paths Options')

    group.add_option('-k', '--k-param', type='int', default=100,\
        help='The number of shortest paths to find. (default=100)')

    group.add_option('','--allow-mult-targets', action='store_true', default=False,\
        help='By default, PathLinker will remove outgoing edges from targets to ensure that there is only one target on each path.  If --allow-mult-targets is specified, these edges are not removed.')

    group.add_option('','--allow-mult-sources', action='store_true', default=False,\
        help='By default, PathLinker will remove incoming edges to sources to ensure that there is only one source on each path.  If --allow-mult-sources is specified, these edges are not removed.')
     
    group.add_option('', '--include-tied-paths', action='store_true', default=False,\
        help='Return the k shortest paths as well as all paths with length equal to that of the kth path.')

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

    # Validate options
    if(opts.PageRank and opts.no_log_transform):
        sys.exit('\nERROR: Options --PageRank and --no-log-transform should not be used together. PageRank weights are probabilities, and must be log-transformed to have an additive interpretation.')

    net = None
    ## Read the network from file
    with open(NETWORK_FILE, 'r') as network_file:
        net = pl.readNetworkFile(network_file, opts.PageRank)

    # Print info about the network
    print(nx.info(net))

    # Operate on only the largest connected component
    if(opts.largest_connected_component):

        conn_comps = nx.weakly_connected_component_subgraphs(net)
        
        # This is the only portion of the program which prevents
        # compatibility between Python 2 & 3. In 2, this object is a
        # generator, but in 3 it is a list. Just check the type and
        # handle accordingly to provide cross-compatibility.
        if(isinstance(conn_comps, types.GeneratorType)):
            net = next(conn_comps)
        elif(isinstance(conn_comps, list)):
            net = conn_comps[0]
        else:
            print("Compatibility error between NetworkX and Python versions. Connected components object from NetworkX does not have acceptable type.")
            exit(-1)

        print("\n Using only the largest weakly connected component:\n"+nx.info(net))

    # Read the sources and targets on which to run PageRank and KSP
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


    # Run PageRank on the network
    # (if opts.PageRank == false, the weights were read from a file above)
    if(opts.PageRank):

        PR_PARAMS = {'q' : opts.q_param,\
                     'eps' : opts.epsilon,\
                     'maxIters' : opts.max_iters}

        print('\nRunning PageRank on net.(q=%f)' %(opts.q_param))

        # The initial weights are entirely on the source nodes, so this 
        # corresponds to a random walk that teleports back to the sources.
        weights = dict.fromkeys(sources, 1.0)
        prFinal = pagerank(net, weights, **PR_PARAMS)
        
        # Write node visitation probabilities
        # (function imported from PageRank)
        writePageRankWeights(prFinal,filename='%s-node-pagerank.txt' % (opts.output))

        # Weight the edges by the flux from the nodes
        pl.calculateFluxEdgeWeights(net, prFinal)

        # Write edge fluxes
        pl.printEdgeFluxes('%s-edge-fluxes.txt' % (opts.output), net)
    
    ## Prepare the network to run KSP

    # Remove improper edges from the sources and targets. If the user runs PageRank 
    # first, then this step will cause the total edge flux in the graph to be less than one.  
    # These transformations are executed by default to prevent them, use the 
    # opts.allow_mult_sources or opts.allow_mult_target arguments.
    if not opts.allow_mult_sources:
        pl.modifyGraphForKSP_removeEdgesToSources(net, sources)
    if not opts.allow_mult_targets:
        pl.modifyGraphForKSP_removeEdgesFromTargets(net, targets)

    # Apply the user specified edge penalty
    pl.applyEdgePenalty(net, opts.edge_penalty)
        
    # Transform the edge weights with a log transformation
    if(not opts.no_log_transform):
        pl.logTransformEdgeWeights(net)

    # Add a super source and super sink. Performed after the
    # transformations so that the edges can be given an additive
    # weight of 0 and thus not affect the resulting path cost.
    pl.modifyGraphForKSP_addSuperSourceSink(
        net, sources, targets, weightForArtificialEdges = 0)
    
    ## Run the pathfinding algorithm
    print('\nComputing the k=%d shortest simple paths.' %(opts.k_param))
    paths = ksp.k_shortest_paths_yen(net, 'source', 'sink', opts.k_param, weight='ksp_weight', clip=not opts.include_tied_paths)

    if len(paths)==0:
        sys.exit('\tERROR: Targets are not reachable from the sources.')

    ## Use the results of KSP to rank edges

    # Un-does the logarithmic transformation on the path lengths to
    # make the path length in terms of the original edge weights
    if not opts.no_log_transform:
        paths = pl.undoLogTransformPathLengths(paths)

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
                edges.append( (t,h,{
                    'ksp_id':k, 
                    'ksp_weight':net.edge[t][h]['ksp_weight'],
                    'path_cost': path[-1][1]}) )

        # Add all new, good edges from this path to the network
        pathgraph.add_edges_from(edges)

        # Each node is ranked by the first time it appears in a path.
        # Identify these by check for any nodes in the graph which do
        # not have 'ksp_id' attribute, meaning they were just added
        # from this path.
        for n in pathgraph.nodes():
            if 'ksp_id' not in pathgraph.node[n]:
                pathgraph.node[n]['ksp_id'] = k


    ## Write out the results to file    

    # Write a list of all edges encountered, ranked by the path they
    # first appeared in.
    kspGraphOutfile = '%sk-%d-ranked-edges.txt' %(opts.output, opts.k_param)
    pl.printKSPGraph(kspGraphOutfile, pathgraph)
    print('\nKSP results are in "%s"' %(kspGraphOutfile))

    # Write a list of all paths found by the ksp algorithm, if
    # requested.
    if(opts.write_paths):
        kspOutfile = '%sk-%d-paths.txt' %(opts.output, opts.k_param)

        pl.printKSPPaths(kspOutfile, paths)
        print('KSP paths are in "%s"' %(kspOutfile))

    print('\nFinished!')


if __name__ == '__main__':
    main(sys.argv)
