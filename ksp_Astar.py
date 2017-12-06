'''
Yen's k shortest loopless paths algorithm, augmented to use A* as a pathfinding
subroutine to improve performance.

For more information, see the following work:
Pathways on Demand: Automated Reconstruction of Human Signaling Networks
Anna Ritz, Christopher L. Poirel, Allison N. Tegge, Nicholas Sharp, Allison Powell, Kelsey Simmons, Shiv D. Kale, and T. M. Murali
Virginia Tech, Blacksburg, VA
Manuscript under review.

Relevant reference:
    Yen, J. Y. (1971). Finding the k shortest loopless paths in a network.
    management Science, 17(11), 712-716.

This code is authored by:
Nicholas Sharp: nsharp3@vt.edu
Anna Ritz: annaritz@vt.edu
Christopher L. Poirel: chris.poirel@gmail.com
T. M. Murali: tmmurali@cs.vt.edu
'''


import sys
from optparse import OptionParser, OptionGroup
import heapq
from collections import defaultdict
from math import isinf

import networkx as nx


# A general implementation of the A* algorithm.
# Compute exact shortest paths in the network, utilizing the fact
# that the heuristic function is monotonic.
# Modeled after the networkx single_source_dijkstra() function
def shortest_path_Astar(net, source, target, heuristicF, weight):

    # Terminology: A node is PROCESSED once we have found the minimum path
    # to it (at that point it gets added to the dist[] dict. A node is SEEN
    # once it it has been encountered as a neighbor to a PROCESSED node.

    # The A* algorithm for monotonic heuristic functions is very similar to
    # Dijkstra's algorithm. The difference is that we key the heap by the
    # distance to the node PLUS the heuristic function, which gives a lower
    # bound on the length of the path the rest of the way to the target.

    # This algorithm follows the convention that if the value of the
    # heuristic is float('inf'), this corresponds to being unreachble.
    # That is, the node is disconnected from the target according the
    # heuristic. This means that float('inf') cannot be used as a valid
    # edge weight, and heuristic functions must conform.

    if source==target:
        return ({source:0}, {source:[source]})

    dist = {}                   # Dictionary of final distances
    preds = {source:None}       # Dictionary of paths
    seen = {source:0}           # Map of seen nodes and their distances
    fringe = []                 # Heap of nodes on the border to process, keyed by heuristic distance
                                # -- Use the format (distance, label)

    heapq.heappush(fringe, (heuristicF(source), (source, 0)))

    # Real-valued edge weights can cause the search to fail due to accumulated
    # error summing along the path. Test with a relative epsilon to catch only
    # the'real' errors.
    REL_EPS = 1e-10

    # Iteratively search the graph outward until we've processed all
    # nodes
    while fringe:

        (heurDist, (currNode, actualDist)) = heapq.heappop(fringe)

        # If we've already processed this node, don't re-process it.
        # This happens because when we see a better path to an already
        # seen node, it's cheaper to leave it in the heap and deal
        # with it here than to remove it (noted below).
        if currNode in dist:
            continue

        # Process this node, this is necessarily the best possible path
        # to it.
        dist[currNode] = actualDist

        # Check for a soluion.
        if currNode == target:
            break

        # Examine all neighbors to this node and consider adding them
        # to the fringe
        currEdges = iter(net[currNode].items())
        for nextNode, edgedata in currEdges:

            # The actual distance to the node from the source (using currNode)
            nextActDist = actualDist + edgedata.get(weight,1)

            # The heuristic function gives a lower bound on the path length
            # to go the rest of the way to the finish from the start
            nextHeurDist = nextActDist + heuristicF(nextNode)

            # If the heuristic returns float('inf'), then the target is
            # necessarily unreachable, so don't expand the search along
            # this edge
            if isinf(nextHeurDist):
                continue

            # If we've already processed the neighbor, then this can't possibly be
            # a better path, assuming the problem is well-formed.
            if nextNode in dist:
                # Verfiy that the graph and heuristic don't break the search property
                if (nextActDist * (1+REL_EPS)) < dist[nextNode]:
                    raise ValueError('Contradictory search path:', 'bad heuristic? negative weights?')

            # If this node hasn't already been processed, we need to consider adding
            # it to the heap. If it's not already in the heap, we must add it. If it
            # is already in the heap, we should only add it if this path to it is an
            # improvement over the previous path. For performance, we leave the old
            # entry in the heap in that case and skip it when it pops
            # out.
            elif nextNode not in seen or nextActDist < seen[nextNode]:
                seen[nextNode] = nextActDist
                heapq.heappush(fringe, (nextHeurDist, (nextNode, nextActDist)))
                preds[nextNode] = currNode

    return (preds, dist)

def build_path_from_preds(preds, source, target):

    pathList = []
    currNode = target
    while currNode != source:
        pathList.append(currNode)
        currNode = preds[currNode]
    pathList.append(currNode)

    pathList.reverse()
    return pathList

# Compute the k shortest acyclic paths in the graph using Yen's
# algorithm. Assumes that this is NOT a multigraph (there is at most one
# edge between any two nodes).
#
# A* is used as the pathfinding subroutine, with the distances in the
# input graph as a heurstic. Because the algorithm computes paths over
# subsets of the initial graph, this heuristic is valid and effective.
#
# Accepts an optional argument 'thresh', which requires that the algorithm
# should not terminate until it has found a path with weight which
# is a factor of 'thresh' larger than that of this initial path. If
# necessary more than k paths are returned to satisfy this threshold.
#
# If clip is set to False, then the method will always return all paths
# of equal length, even if doing so returns more than k paths. That is,
# if k = 10 and the 10th and 11th paths both have the same weight, the
# 11th path will be included.
#
# If the graph contains n < k paths, n paths will be returned.
def k_shortest_paths_yen(G, source, target, k=1, weight='weight', thresh=None, clip=True):
    net = G.copy()

    if source==target:
        return [(source, 0.0)]
    if source not in net or target not in net:
        return []

    # Compute the original distance from the source to use for the heuristic
    # function.
    net.reverse(copy=False)
    minDists = nx.single_source_dijkstra_path_length(net, target, weight=weight)
    net.reverse(copy=False)

    # Heuristic function for the A* search, using distances
    # in the original network. Note that this function is monotone.
    def heuristicF(u):
        if u in minDists:
            return minDists[u]
        else:
            return float("inf")    # corresponds to disconnected


    # Compute the initial shortest path to initialize Yen's
    shortestSourceDists, shortestSourcePaths = nx.single_source_dijkstra(net, source, target, weight=weight)

    if target not in shortestSourcePaths:
        return []

    # prevPath stores the most recently computed shortest path from source to target.
    # it is a list of pairs, where each pair contains two elements:
    #    (curr node in path, cost of path up to current node)
    prevPath = [(x, shortestSourceDists[x]) for x in shortestSourcePaths[target]]

    # stores k shortest paths in the order they are discovered
    shortestPaths = [prevPath]
    candidates = []

    # A cache mapping prefixes of accepted paths to the next node after
    # the prefix. Used to avoid scanning all previous paths many time,
    # which otherwise dominates runtime.
    prefixCache = defaultdict(list)
    for i in range(1,len(prevPath)):
        prefixCache[tuple(prevPath[:i])].append(prevPath[i][0])

    if thresh == None:
        threshSatisfied = True
    else:
        threshSatisfied = False
        initialPathCost = prevPath[-1][1]

    # Find the paths. This loop is structured in a slightly unusual way
    # to support the 'clip=False' option, which requires looking ahead
    # to the k+1 path.
    currk = 2
    while True:

        # Check for completion. The loop actually terminates below after
        # checking 'clip' if needed.
        readyToTerminate = currk > k and threshSatisfied
        if(readyToTerminate and clip):
            break;

        # Save the edges that are hidden so we can add then back in
        # after the iteration.
        hiddenEdges = []

        # Process each node of the most recently found path, computing
        # the shortest path that deviates at that node and adding it
        # to the candidate heap.
        for i in range(len(prevPath)-1):
            x, xcost = prevPath[i]

            # hide edges incoming to x until iteration currk is over to
            # avoid finding cycles. Note that this effect is cumlative,
            # meaning that while processing the current node in the path,
            # all incoming edges to this node and all previous nodes
            # have been hidden.
            for u,v, edata in net.in_edges(x, data=True):
                hiddenEdges.append( (u,v, edata) )
                net.remove_edge(u,v)

            # for each previously-found shortest path P_j with the same first i nodes as
            # the first i nodes of prevPath, hide the edge from x to the i+1 node in P_j
            # to ensure we don't re-find a previously found path.
            # Lookup the prefixes in a cache to disallow them. Requires more memory 
            # to store the cache, but saves scanning the list of found paths
            for repNode in prefixCache[tuple(prevPath[:i+1])]:
                if net.has_edge(x, repNode):
                    hiddenEdges.append( (x, repNode, net.get_edge_data(x,repNode)) )
                    net.remove_edge(x,repNode)


            # concatenate prevPath[:i+1] and shortest path from x->target,
            # and add this path to candidates

            # Find the shortest path using A*.
            (spurPreds, spurDists) = shortest_path_Astar(net, x, target, heuristicF, weight)

            # Short circuit if the target node was unreachable, which is expected
            # to happen as we remove edges.
            if target not in spurDists:
                continue
            spurPath = build_path_from_preds(spurPreds, x, target)

            # Add this candidate path to the heap if its not already there. The distances
            # need to be incremented with the cost of the prefix
            spurDist = spurDists[target] + xcost
            spurPath = [(n, spurDists[n] + xcost) for n in spurPath]
            newCandidate = (spurDist, prevPath[:i] + spurPath)
            if newCandidate not in candidates:
                r = heapq.heappush(candidates, newCandidate)


        # Only need to store k-currk+1 shortest paths in the heap. This
        # gets expensive so don't do it every iteration.
        # We can only safely do this if we don't have a threshold
        # to satisfy.
        if threshSatisfied and len(candidates)>(k-currk+1) and (currk % 100) == 0:
            keepCandidates = heapq.nsmallest(k-currk+2, candidates, key=lambda x: x[0])
            candidates = keepCandidates
            heapq.heapify(candidates)

        # terminate early if there are no more paths from source to vertex
        if len(candidates)==0:
            break

        # Accept the shortest path on the candidates heap, which is
        # necessarily the next shortest path.
        newShortest = heapq.heappop(candidates)[1]

        # Check for termination if clip=False (if readyToTerminate=True
        # and clip=True we would have already exited at the check above)
        if readyToTerminate:
            # If this path is longer than the previous one than we can
            # exit even if clip=True
            if newShortest[-1][1] != prevPath[-1][1]:
                break

        # Add this to the list of prefixes for efficient lookup later
        for i in range(1,len(newShortest)):
            if not newShortest[i][0] in prefixCache[tuple(newShortest[:i])]:
                prefixCache[tuple(newShortest[:i])].append(newShortest[i][0])

        shortestPaths.append(newShortest)
        prevPath = newShortest

        # add the hidden edges back to the network before the next iteration
        net.add_edges_from(hiddenEdges)

        # Check if this new shortest path causes the threshold to be satisfied
        if ((not threshSatisfied) and newShortest[-1][1] >= initialPathCost*thresh):
            threshSatisfied = True

        currk += 1

    return shortestPaths

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
        outf.write('%d\t%0.5e\t%s\n' %(k, length, '|'.join(pathNodes) ))
    outf.close()
    return

# Main method, so this can be used on the command line
def main(args):

    usage = '''
ksp_Astar.py [options] NETWORK SOURCE TARGET
REQUIRED arguments:
    NETWORK - A tab-delimited file with one directed interaction per line. Each
        line should have at least 2 columns: tail, head. Edges are directed from
        tail->head. This file can optionally have a third column specifying the
        edge weight

    SOURCE - A node name used as the source node for the pathfinding

    TARGET - A node name used as the target node for the pathfinding

'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('-o', '--output', type='string', default='paths.txt', metavar='STR',\
        help='Filename to print the resulting weights. (default="paths.txt")')

    parser.add_option('-k', '--k-param', action='store', type='int', default=200,\
        help='The number of shortest paths to find (default=200)')


    # Parse the command line arguments
    (opts, args) = parser.parse_args()

    # Get the required arguments
    num_req_args = 3
    if len(args)!=num_req_args:
        parser.print_help()
        sys.exit('\nERROR: ksp_Astar.py requires %d positional arguments, %d given.' %(num_req_args, len(args)))
    NETWORK_FILE = args[0]
    SOURCE = args[1]
    TARGET = args[2]


    ## Read in the graph from file
    net = nx.DiGraph()

    # Read the network file
    print('\nReading the network from %s\n' %(NETWORK_FILE))
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

        # If no weight is given for the edge, assign it a weight of 1.
        eWeight = 1
        if(len(items) > 2):
            eWeight = float(items[2])

        net.add_edge(id1, id2, weight=eWeight)

    if not (SOURCE in net and TARGET in net):
        print("ERROr: SOURCE and TARGET must both be nodes in the network")
        exit()

    ## Run ksp
    print("Computing k = %d shortest paths\n"%(opts.k_param))
    paths = k_shortest_paths_yen(net, SOURCE, TARGET, k=opts.k_param, weight='weight')

    ## Print the result
    print("Writing results to " + opts.output + "\n")
    printKSPPaths(opts.output, paths)

if __name__=='__main__':
    main(sys.argv)
