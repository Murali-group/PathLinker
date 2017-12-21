def write_network_file(G, output_file_handle):
    output_file_handle.write("# Tail Node\tHead Node\tWeight\n")
    for edge in G.edges(data=True): 
        line = \
            str(edge[0]) + "\t" + \
            str(edge[1]) + "\t" + \
            str(edge[2]['weight'] + "\n")

        output_file_handle.write(line)


def get_source_set(nodes_file_handle):
    sources = set()
    for line in nodes_file_handle:
        if not is_comment_line(line): 
            tokens = tokenize(line)
            if tokens[1] in ['source', 'receptor']:
                sources.add(tokens[0])
    return sources


def get_target_set(nodes_file_handle):
    targets = set()
    for line in nodes_file_handle:
        if not is_comment_line(line): 
            tokens = tokenize(line)
            if tokens[1] in ['target', 'tr', 'tf']:
                targets.add(tokens[0])
    return targets


def get_node_set(network_file_handle):
    nodes = set()
    for line in network_file_handle:
        if not is_comment_line(line): 
            tokens = tokenize(line)
            nodes.add(tokens[0])
            nodes.add(tokens[1])
    return nodes


def get_edge_set(network_file_handle):
    edges = set()
    for line in network_file_handle:
        if not is_comment_line(line): 
            tokens = tokenize(line)
            edge = get_edge(tokens)
            edges.add(edge)
    return edges


def get_ranked_nodes(ranked_edges):
    """
    Given: a list of sets of edges, where all nodes in the set have
    the same rank

    Outputs: A list of sets of nodes, where all nodes in the set have
    the same rank
    """
    ranked_nodes = []
    overall_set = set()
    current_set = set()

    for edge_set in ranked_edges:
        for edge in edge_set:
            for node in edge:
                if node not in overall_set:
                    overall_set.add(node)
                    current_set.add(node)

        ranked_nodes.append(current_set)
        current_set = set()

    return ranked_nodes
            

def parse_ranked_edges(ranked_edges_handle):
    ranked_edges = []
    current_set = set()
    prev_rank = 1

    for line in ranked_edges_handle:
        if not is_comment_line(line): 
            tokens = tokenize(line)
            edge = get_edge(tokens)
            rank = int(get_edge_rank(tokens))

            if rank == prev_rank:
                current_set.add(edge)
            else:
                ranked_edges.append(current_set)
                current_set = set()
                current_set.add(edge)
                prev_rank = rank

    ranked_edges.append(current_set)

    return ranked_edges


def parse_ranked_paths(ranked_paths_handle):
    ranked_paths = []
    prev_rank = 1

    for line in ranked_paths_handle:
        if not is_comment_line(line):
            tokens = tokenize(line)
            path = get_path(tokens)
            ranked_paths.append(path)

    return ranked_paths


def is_comment_line(line):
    if line.lstrip().startswith("#"):
        return True
    return False


def tokenize(line):
    return line.split("\t")


def get_edge(tokens):
    return (tokens[0], tokens[1])


def get_edge_rank(tokens):
    return tokens[2]


def get_path(tokens):
    return (tokens[2].strip(), tokens[1])


def get_path_rank(tokens):
    return tokens[1]
