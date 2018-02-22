'''
An implementation of the standard pagerank algorithm, using iterative convergence.

Citation for this work:
    Poirel, C. L., Rodrigues, R. R., Chen, K. C., Tyson, J. J., & Murali, T. M. 
    (2013). Top-down network analysis to drive bottom-up modeling of physiological 
    processes. Journal of Computational Biology, 20(5), 409-418.

Relevant reference:
    Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank citation 
    ranking: Bringing order to the web.

This code is authored by:
Nicholas Sharp: nsharp3@vt.edu
Anna Ritz: annaritz@vt.edu
Christopher L. Poirel: chris.poirel@gmail.com
T. M. Murali: tmmurali@cs.vt.edu
"""
'''

# Imports
import networkx as nx
import sys
from optparse import OptionParser, OptionGroup

'''
Run the PageRank algorithm.

Inputs:
    + net           - A NetworkX graph to run the pagerank algorithm on
    + weights       - (Optional) Initial weights for the graph. These are used both
                      to initialize the algorithm and to provide destination probatilities
                      during the teleportation step. If not given, uniform weights are used.
    + q             - The teleportation probability for the PageRank algorithm
                      (default = 0.5)
    + eps           - A RELATIVE convergence threshold for iteration (default = 0.01)
    + maxIters      - The maximum number of iterations to perform (default = 500)
    + verbose       - Print extra information (default = False)
    + weightName    - What property key holds the weights in the graph (default = weight)

Outputs:
    + currVisitProb - A dictionary of the final node probabilities at the end of the
      iteration process.

'''
def pagerank(net, weights={}, q=0.5, eps=0.01, maxIters=500, verbose=False, weightName = 'weight'):

    incomingTeleProb = {}   # The node weights when the algorithm begins, also used as 
                            # teleport-to probabilities

    prevVisitProb = {}      # The visitation probability in the previous iteration

    currVisitProb = {}      # The visitation probability in the current iteration

    N = net.number_of_nodes()

    # Find the incoming teleportation probability for each node, which
    # is also used as the initial probabilities in the graph. If no
    # node weights are passed in, use a uniform distribution.

    totWeight = sum([w for v,w in weights.items()])

    # If no incoming weights are given, use a uniform weighting.
    if totWeight==0:
        incomingTeleProb = dict.fromkeys(net, 1.0/N)
        prevVisitProb = incomingTeleProb.copy()
        currVisitProb = incomingTeleProb.copy()

    # If weights are given, apply two transformations
    #   - Add a small incoming teleportation probability to every node to
    #     ensure that the graph is strongly connected
    #   - Normalize the weights to sum to one: these are now probabilities.
    else:

        # Find the smallests non-zero weight in the graph. 
        minPosWeight = 1.0
        for v, weight in weights.items():
            if weight==0:
                continue
            minPosWeight = min(minPosWeight, 1.0*weight/totWeight)

        # The epsilon used as the added incoming teleportation
        # probabilty is 1/1000000 the size of the smallest weight given
        # so that it does not impact results.
        smallWeight = minPosWeight/(10**6)

        # Explicitly calculate the weight, including the added
        # teleportation probabilitiy.
        for v in net.nodes_iter():
            weight = weights.get(v, 0.0) # return weights[v], 0 otherwise
            incomingTeleProb[v] = 1.0*(weight + smallWeight)/(totWeight + smallWeight*N)
        prevVisitProb = incomingTeleProb.copy()
        currVisitProb = incomingTeleProb.copy()
    # END if/else that initializes teleportation probabilities.

    # Cache out-degree of all nodes to speed up the update loop
    outDeg = {}
    zeroDegNodes = set()
    for v in net.nodes_iter():
        outDeg[v] = 1.0*net.out_degree(v, weight= weightName)
        if outDeg[v]==0:
            zeroDegNodes.add(v)

    # Apply a standard iterative convergence procedure to find the
    # pagerank node weights. See Page et. al. reference for a general
    # explantion.
    iters = 0
    finished = False
    while (not finished):

        iters += 1
        prevVisitProb = currVisitProb.copy()
        maxDiff = 0

        # In the basic formulation, nodes with degree zero ("dangling
        # nodes") have no weight to send a random walker if it does not
        # teleport. We consider a walker on one of these nodes to
        # teleport with uniform probability to any node. Here we compute
        # the probability that a walker will teleport to each node by
        # this process.
        zSum = sum([prevVisitProb[x] for x in zeroDegNodes]) / N

        # Calculate a new visitation probability for each node
        for v in net.nodes_iter():

            # Calculate the total probability of a walker entering this
            # node from any of the neighbors
            eSum = 0
            for u in net.predecessors_iter(v):
                w_uv = 1.0*net[u][v][weightName]
                eSum += w_uv/outDeg[u] * prevVisitProb[u]

            # The total probability of a walker entering this node is
            # the sum of three components
            #   - Walking in from an edge
            #   - Teleporting in
            #   - Teleporting in (from a dangling node)

            # The second term here represents the latter two terms,
            # evaluated analytically to avoid creating an excessive
            # number of edges in the graph.
            currVisitProb[v] = q*incomingTeleProb[v] + (1-q)*(eSum + zSum)

            # Keep track of the maximum RELATIVE difference between this
            # iteration and the previous to test for convergence
            maxDiff = max(maxDiff, abs((prevVisitProb[v] - currVisitProb[v]) / 
                                         currVisitProb[v]))

        # Print statistics on the iteration
        if verbose:
            print('\tIteration %d, max difference %f' %(iters, maxDiff))
            if maxDiff < eps:
                print('PageRank converged after %d iterations, max difference %f.' %(iters, maxDiff))

        # Give a warning if termination happens by the iteration cap,
        # which generally should not be expcted.
        if iters >= maxIters:
            print('WARNING:PageRank terminated because max iters (%d) was reached.' %(maxIters))

        # Test for termination, either due to convergence or exceeding
        # the iteration cap
        finished = (maxDiff<eps) or (iters>=maxIters)

    return currVisitProb


def writePageRankWeights(final, filename=None):
    '''
    Write the resulting weights from a pagerank run to file for later use
    '''
    ostream = sys.stdout
    if filename!=None:
        ostream = open(filename, 'w')
    ostream.write('#nodeID\tfinal\n')
    for n, w in sorted(final.items(), key=lambda x: x[1], reverse=True):
        ostream.write('%s\t%s\n' %(n, final[n]))
    if ostream!=sys.stdout:
        ostream.close()

def main(args):
    '''
    Main method, so this can be used on the command line
    '''

    usage = '''
PageRank.py [options] NETWORK
REQUIRED arguments:
    NETWORK - A tab-delimited file with one directed interaction per line. Each
        line should have at least 2 columns: tail, head. Edges are directed from
        tail->head. This file can optionally have a third column specifying the
        edge weight

'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('-o', '--output', type='string', default='pageranks.txt', metavar='STR',\
        help='Filename to print the resulting weights. (default="pageranks.txt")')

    parser.add_option('-v', '--verbose', action='store_true', default=False,\
        help='Print statistics about each iteration and other information')

    # Random Walk Parameter Group
    group = OptionGroup(parser, 'Random Walk Options')

    group.add_option('-q', '--q-param', action='store', type='float', default=0.5,\
        help='The value of q indicates the probability that the random walker teleports back to a source node during the random walk process. (default=0.5)')

    group.add_option('-e', '--epsilon', action='store', type='float', default=0.01,\
        help='A small value used to test for relative convergence of the iteration. (default=0.01)')

    group.add_option('', '--max-iters', action='store', type='int', default=500,\
        help='Maximum number of iterations to run the PageRank algorithm. (default=500)')

    parser.add_option('', '--tele-weights', type='string', default=None, metavar='STR',\
        help='Incoming teleportation weights for each node, in a tab-separated file ' + \
        'with the form "nodename[tab]weight". If not given, uniform weights are used.')
    
    parser.add_option_group(group)

    # Parse the command line arguments
    (opts, args) = parser.parse_args()

    # Get the required arguments
    num_req_args = 1
    if len(args)!=num_req_args:
        parser.print_help()
        sys.exit('\nERROR: PageRank.py requires %d positional arguments, %d given.' %(num_req_args, len(args)))
    NETWORK_FILE = args[0]


    ## Read in the graph from file
    net = nx.DiGraph()

    # Read the network file
    print('\nReading the network from %s' %(NETWORK_FILE))
    infile = open(NETWORK_FILE, 'r')
    for line in infile:
        items = [x.strip() for x in line.rstrip().split('\t')]

        # Skip empty lines or those beginning with '#' comments
        if (line == '') or (line[0] =='#'):
            continue

        id1 = items[0]
        id2 = items[1]

        # If no weight is given for the edge, assign it a weight of 1.
        eWeight = 1
        if(len(items) > 2):
            eWeight = float(items[2])

        net.add_edge(id1, id2, weight=eWeight)


    ## Read teleportation weights if given
    teleProbs = {} # Node --> weight
    if opts.tele_weights != None:
        
        print('\nReading incoming teleportation probabilities from %s' %(opts.tele_weights))
        infile = open(opts.tele_weights, 'r')
        
        for line in infile:
            items = [x.strip() for x in line.rstrip().split('\t')]

            # Skip empty lines or those beginning with '#' comments
            if (line == '') or (line[0] =='#'):
                continue

            node = items[0]
            weight = float(items[1])

            if not node in net:
                print("Error: Node %s from teleportation probability file is not in graph."%(node))
                exit(-2)

            teleProbs[node] = weight
    
    ## Run PageRank
    finalProbs = pagerank(net, weights=teleProbs,
            q=opts.q_param, eps=opts.epsilon, maxIters=opts.max_iters, verbose=opts.verbose)


    ## Print the result
    print("\nWriting results to " + opts.output)
    writePageRankWeights(finalProbs, filename=opts.output)


if __name__=='__main__':
    main(sys.argv)
