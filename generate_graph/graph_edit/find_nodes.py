#!/usr/bin/env python3

import networkx
from graph_cleaning_exception import GraphCleaningDFSException


def find_nodes_at_distance(ontology, dst, nodes, graph):
    """
    Find nodes at a specific distance from the root node in the graph.

    Args:
        ontology (str): Ontology identifier for the layer attribute.
        dst (int): Target distance from the root node.
        nodes (iterable): Nodes to search within.
        graph (networkx.DiGraph): Directed graph.

    Returns:
        list: Nodes at the specified distance from the root node.
    """
    result_nodes = [
        node for node in nodes
        if node in graph.nodes and f"Layer_of{ontology}" in graph.nodes[node]
        and graph.nodes[node][f"Layer_of{ontology}"] == dst
    ]
    return result_nodes

def find_farthest_node(ontology, start_node, graph):
    """
    Find the farthest node in a directed graph using DFS.

    Args:
        ontology (str): Ontology identifier for the layer attribute.
        start_node: Starting node for the DFS traversal.
        graph (networkx.DiGraph): Directed graph.

    Returns:
        The node with the greatest distance in the ontology's layer attribute.

    Raises:
        GraphCleaningDFSException: If a node does not have the required layer attribute.
    """
    farthest_node = start_node
    for edge in networkx.edge_dfs(graph, source=start_node):
        target_node = edge[1]
        layer_key = f"Layer_of{ontology}"

        if layer_key not in graph.nodes[target_node]:
            raise GraphCleaningDFSException(start_node, ontology)

        if graph.nodes[target_node][layer_key] > graph.nodes[farthest_node][layer_key]:
            farthest_node = target_node

    return farthest_node


def mark_lowest_nodes(pivots, graph):
    """
    Mark the lowest nodes in the graph by traversing from nodes with 'UniprotID' attributes.

    Args:
        pivots (list): List of ontology identifiers.
        graph (networkx.DiGraph): Directed graph.

    Returns:
        networkx.DiGraph: Updated graph with marked lowest nodes.
    """
    for ontology in pivots:
        layer_key = f"Layer_of{ontology}"
        for node in graph.nodes():
            if layer_key in graph.nodes[node] and 'UniprotID' in graph.nodes[node] and graph.nodes[node]['UniprotID'] != "pruchozi":
                farthest_node = find_farthest_node(ontology, node, graph)

                if "UniprotID" not in graph.nodes[farthest_node]:
                    graph.nodes[farthest_node]["UniprotID"] = "lowest"
    return graph

def print_graph_info(ontology, distance, nodes, graph, result_nodes):
    print(f"""
Ontology: {ontology}, Distance: {distance}
Nodes checked: {nodes}
Nodes in graph: {list(graph.nodes)}
Found nodes at distance: {result_nodes}
""")



if __name__ == "__main__":

    G = networkx.DiGraph()
    G.add_nodes_from([(2, {'Layer_ofGO:0003674': 1}), (22, {'Layer_ofGO:0003674': 1}), ('GO:0003674', {'Layer_ofGO:0003674': 0}), ('UniprotID', {'UniprotID': 'P05067', 'Layer_ofGO:0003674': 2}), (5, {'Layer_ofGO:0003674': 4}), (6, {'Layer_ofGO:0003674': 3}), (7, {'Layer_ofGO:0003674': 5, 'UniprotID': 'pruchozi'}), (3, {'Layer_ofGO:0003674': 2}), (8, {'Layer_ofGO:0003674': 4})])
    G.add_edges_from([(2, 'GO:0003674'), ('UniprotID', 2), ('UniprotID', 5), (5, 6), (5, 7), (6, 3), (7, 8), (3, 2), (8, 6)])

    ontology = 'GO:0003674'
    distance = 1

    nodes_at_distance_from_root = find_nodes_at_distance(ontology, distance, G.nodes(), G)

    print_graph_info(ontology, distance, G.nodes(), G, nodes_at_distance_from_root)
