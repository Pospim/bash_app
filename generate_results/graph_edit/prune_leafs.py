#!/usr/bin/env python3

import networkx
from find_nodes import find_nodes_at_distance

def to_remove(removable, neighbors, leaf):
    removable_set = set(removable+[leaf])
    return set(neighbors).issubset(removable_set)

def update_dst(update_flag, neighbor, dst, ontology, graph):
    """
    Update the distance and modify flag based on the neighbor's layer.

    Args:
        update_flag (bool): Current flag indicating distance modification.
        neighbor: Neighbor node being checked.
        dst (int): Current distance value.
        ontology (str): Ontology identifier.
        graph (networkx.DiGraph): Directed graph.

    Returns:
        tuple: Updated flag and distance value.
    """
    layer_key = f'Layer_of{ontology}'
    if graph.nodes[neighbor][layer_key] >= dst:
        return True, graph.nodes[neighbor][layer_key]
    return update_flag, dst

def process_leaves(removable, leaves, ontology, dst, graph, reverse_graph):
    """
    Process removable nodes and update leaves and distance.

    Args:
        removable (list): Nodes to be removed.
        leaves (list): Current list of leaves.
        ontology (str): Ontology identifier.
        dst (int): Current distance value.
        graph (networkx.DiGraph): Directed graph.
        reverse_graph (networkx.DiGraph): Reverse of the directed graph.

    Returns:
        tuple: Updated leaves and distance.
    """
    updated = False
    while removable:
        leaf = removable.pop()
        for neighbor in graph.neighbors(leaf):
            if graph.in_degree(neighbor) > 1:
                neighbors = list(reverse_graph.neighbors(neighbor))
                if neighbor not in leaves and to_remove(removable, neighbors, leaf):
                    leaves.append(neighbor)
                    updated, dst = update_dst(updated, neighbor, dst, ontology, graph)
            else:
                leaves.append(neighbor)
                updated, dst = update_dst(updated, neighbor, dst, ontology, graph)
        graph.remove_node(leaf)
        reverse_graph.remove_node(leaf)
        if leaf in leaves:
            leaves.remove(leaf)

    if not updated:
        dst -= 1
    return leaves, dst

def clean_graph_layer(ontology, leaves, max_dst, graph, reverse_graph):
    """
    Clean the graph by removing layers without 'UniprotID'.

    Args:
        ontology (str): Ontology identifier.
        leaves (list): Current list of leaves.
        max_dst (dict): Maximum distance for each ontology.
        graph (networkx.DiGraph): Directed graph.
        reverse_graph (networkx.DiGraph): Reverse of the directed graph.

    Returns:
        tuple: Cleaned graph and updated leaves.
    """
    dst = max_dst[ontology]
    removable = []
    while dst > 0:

        leaves_in_distance = find_nodes_at_distance(ontology, dst, leaves, graph)
        for leaf in leaves_in_distance:
            if 'UniprotID' in graph.nodes[leaf].keys():
                pass
            elif leaf not in removable:
                removable.append(leaf)
        leaves, dst = process_leaves(removable, leaves, ontology, dst, graph, reverse_graph)

    return graph, leaves


def clean_all_layers(graph_dict, leaves, max_dst, ontologies):
    """
    Clean graphs for all ontologies by removing unnecessary layers.

    Args:
        graph_dict (dict): Dictionary of ontology-specific graphs.
        leaves (dict): Leaf nodes for each ontology.
        max_dst (dict): Maximum distance for each ontology.
        ontologies (list): List of ontology identifiers.

    Returns:
        tuple: Updated graph dictionary and leaves.
    """
    for ontology in ontologies:
        reverse_graph = graph_dict[ontology].reverse()
        graph_dict[ontology], leaves[ontology] = clean_graph_layer(
            ontology, leaves[ontology], max_dst, graph_dict[ontology], reverse_graph
            )
    return graph_dict, leaves


def print_graph_info(graph_dict, leaves, ontologies):
    for ontology in ontologies:
        graph = graph_dict[ontology]
        if 'name' in graph.nodes[ontology]:
            print(f"Ontology: {ontology} - {graph.nodes[ontology]['name']}")
        else:
            print(f"Ontology: {ontology}")
        print(f"Leaves in graph: {leaves[ontology]}")
        print(f"Nodes in cleaned graph: {list(graph.nodes)}")

if __name__ == "__main__":

    pivot = ['GO:0003674'] #molekular function
    """
             'GO:0008150', #bioloical process
             'GO:0005575'] #celular component
    """
    G = networkx.DiGraph()
    G.add_nodes_from([(2, {'Layer_ofGO:0003674': 1}), ('GO:0003674', {'Layer_ofGO:0003674': 0}), ('UniprotID', {'UniprotID': 'P05067', 'Layer_ofGO:0003674': 2}), (5, {'Layer_ofGO:0003674': 4}), (6, {'Layer_ofGO:0003674': 3}), (7, {'Layer_ofGO:0003674': 5, 'UniprotID': 'pruchozi'}), (3, {'Layer_ofGO:0003674': 2}), (8, {'Layer_ofGO:0003674': 4})])
    G.add_edges_from([(2, 'GO:0003674'), ('UniprotID', 2), ('UniprotID', 5), (5, 6), (5, 7), (6, 3), (7, 8), (3, 2), (8, 6)])

    G.add_nodes_from([(9, {'Layer_ofGO:0003674': 6}), (10, {'Layer_ofGO:0003674': 6})])
    G.add_edges_from([(9, 7), (10, 7)])

    #networkx.draw(G, with_labels=True)
    #plt.show()

    clean_Graph = {'GO:0003674': G }


    leaves = {'GO:0003674': ['UniprotID', 10, 9]}
    most_distance = {'GO:0003674': 6}

    print("Before calculation: ", G.nodes())
    clean_Graph, leaves = clean_all_layers(clean_Graph, leaves, most_distance, pivot)
    print_graph_info(clean_Graph, leaves, pivot)
