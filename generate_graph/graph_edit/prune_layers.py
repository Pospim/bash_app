#!/usr/bin/env python3

import networkx
from find_nodes import find_nodes_at_distance
from concurrent.futures import ThreadPoolExecutor

def has_uniprot_id(node, graph):
    return 'UniprotID' in graph.nodes[node]

def uniprot_id_in_layer(nodes, graph):
    """
    Check if any node in a layer contains a 'UniprotID'.

    Args:
        nodes (list): Nodes in the layer.
        graph (networkx.DiGraph): Graph containing the nodes.

    Returns:
        bool: True if any node in the layer contains a 'UniprotID', False otherwise.
    """
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(lambda node: has_uniprot_id(node, graph), nodes))
    return any(results)


def no_lower_layer_edge(node, ontology, graph):
    """
    Check if a node has no edges originating from a lower layer.

    Args:
        node: Node to check.
        ontology (str): Ontology identifier for the layer attribute.
        graph (networkx.DiGraph): Graph containing the node.

    Returns:
        bool: True if no edges originate from a lower layer, False otherwise.
    """
    for neighbor in graph.reverse().neighbors(node):
        if graph.nodes[node][f'Layer_of{ontology}'] >= graph.nodes[neighbor][f'Layer_of{ontology}']:
            return False
    return True

def add_leaves(leaves, node, ontology, graph):
    """
    Add new leaves to the list of leaves if they meet certain conditions.

    Args:
        leaves (list): Current list of leaves.
        node: Node being checked.
        ontology (str): Ontology identifier.
        graph (networkx.DiGraph): Graph containing the node.

    Returns:
        list: Updated list of leaves.
    """
    for neighbor in graph.neighbors(node):
        if neighbor not in leaves and no_lower_layer_edge(neighbor, ontology, graph):
            leaves.append(neighbor)
    for neighbor in graph.reverse().neighbors(node):
        if neighbor not in leaves and no_lower_layer_edge(neighbor, ontology, graph):
            leaves.append(neighbor)
    return leaves

def prune_layers_without_uniprot(ontology, graph, max_dst, leaves):
    """
    Remove layers that do not contain a 'UniprotID' and update the graph and leaf nodes.

    Args:
        ontology (str): Ontology identifier for the layer attribute.
        graph (networkx.DiGraph): Graph to modify.
        max_dst (dict): Maximum distance for each ontology.
        leaves (list): Current list of leaf nodes.

    Returns:
        tuple: Updated graph, leaves, and maximum distance for the ontology.
    """
    dst = max_dst[ontology]
    while dst >= 0:
        nodes = find_nodes_at_distance(ontology, dst, graph.nodes(), graph)

        if uniprot_id_in_layer(nodes, graph):
            max_dst[ontology] = dst
            break
        else:
            for node in nodes:
                leaves = add_leaves(leaves, node, ontology, graph)
                if node in leaves:
                    leaves.remove(node)
                graph.remove_node(node)
            dst -= 1
    return graph, leaves, max_dst[ontology]

def get_clean_graph(leaves, max_dst, subgraph_ontology, ontologies, graph):
    """
    Create clean graphs for each ontology by pruning unnecessary layers.

    Args:
        leaves (dict): Leaf nodes for each ontology.
        max_dst(dict): Maximum distance for each ontology.
        subgraph_ontology (dict): Subgraphs for each ontology.
        ontologies (list): List of ontology root nodes.
        graph (networkx.DiGraph): Original graph.

    Returns:
        tuple: Cleaned graphs, updated leaves, and updated max distances.
    """
    clean_graphs = {}

    for ontology in ontologies:

        frozen_subgraph = graph.subgraph(subgraph_ontology[ontology])
        clean_graph = networkx.MultiDiGraph()
        clean_graph.add_edges_from(frozen_subgraph.edges())

        for node in clean_graph:
            clean_graph.nodes[node].update(graph.nodes[node])

        clean_graphs[ontology], leaves[ontology], max_dst[ontology] = prune_layers_without_uniprot(
            ontology, clean_graph, max_dst, leaves[ontology]
        )

    return clean_graphs, leaves, max_dst

def connect_pivot(pivots, graph):
    """
    Add a new node to the graph and connect it to all nodes in the pivot list.

    Args:
        pivots (list): List of pivot nodes to be connected.
        graph (networkx.DiGraph): Directed graph to modify.

    Returns:
        tuple: Updated graph and the new pivot list containing the connecting node.
    """
    # Add a new node with a unique attribute
    graph.add_node("new_connecting_node", type="temporary_node", name="Temporary connecting node")

    # Add edges from each pivot to the new node
    for pivot in pivots:
        graph.add_edge(pivot, "new_connecting_node")

    # Update pivot list to include only the new node
    pivots = ["new_connecting_node"]
    return graph, pivots

def print_graph_info(clean_graphs, leaves, max_distance, ontologies, graph):
    """
    Print information about the cleaned graphs and their properties.

    Args:
        clean_graphs (dict): Cleaned graphs for each ontology.
        leaves (dict): Leaf nodes for each ontology.
        max_distance (dict): Maximum distance for each ontology.
        ontologies (list): List of ontology root nodes.
        graph (networkx.DiGraph): Original graph.
    """
    for ontology in ontologies:
        print(graph.nodes())
        if 'name' in graph.nodes[ontology]:
            print(f"Ontology: {ontology} - {graph.nodes[ontology]['name']}")
        else:
            print(f"Ontology: {ontology}")
        print(f"Nodes in ontology graph: {clean_graphs[ontology].nodes()}")
        print(f"Leaves in ontology graph: {leaves[ontology]}")
        print(f"Max distance in ontology graph: {max_distance[ontology]}")


if __name__ == "__main__":
    pivot = ['GO:0003674'] #molekular function
    """
             'GO:0008150', #bioloical process
             'GO:0005575'] #celular component
    """

    G = networkx.DiGraph()
    G.add_nodes_from([(2, {'Layer_ofGO:0003674': 1}), ('GO:0003674', {'Layer_ofGO:0003674': 0}), ('UniprotID', {'UniprotID': 'P05067', 'Layer_ofGO:0003674': 2}), (5, {'Layer_ofGO:0003674': 4}), (3, {'Layer_ofGO:0003674': 2}), (6, {'Layer_ofGO:0003674': 3}), (7, {'Layer_ofGO:0003674': 5, 'UniprotID': 'pruchozi'}), (8, {'Layer_ofGO:0003674': 4}),
                      (9, {'Layer_ofGO:0003674': 6}), (10, {'Layer_ofGO:0003674': 7}), (11, {'Layer_ofGO:0003674': 7}), (12, {'Layer_ofGO:0003674': 8})])
    G.add_edges_from([(2, 'GO:0003674'), ('UniprotID', 2), ('UniprotID', 5), (5, 6), (5, 7), (3, 2), (6, 3), (7, 8), (8, 6), (9,7), (10,9), (11,9), (12,11), (12,10)])

    print(G.nodes())


    subgraph_ontology = {'GO:0003674': G.nodes()}
    leaves = {'GO:0003674': ['UniprotID', 12]}
    most_distance = {'GO:0003674': 8}

    #networkx.draw(G, with_labels=True)
    #plt.show()
    clean_Graph, leaves, most_distance = get_clean_graph(leaves, most_distance, subgraph_ontology, pivot, G)
    print_graph_info(clean_Graph, leaves, most_distance, pivot, G)
