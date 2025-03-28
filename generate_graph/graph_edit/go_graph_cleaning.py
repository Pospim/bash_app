#!/usr/bin/env python3

import obonet
from add_uniprot_to_node import map_uniprot_to_graph
from add_layers_to_node import get_leaves_by_ontology
from find_nodes import mark_lowest_nodes
from prune_layers import get_clean_graph
from prune_leafs import clean_all_layers
from graph_cleaning_exception import GraphCleaningDFSException
from prune_layers import connect_pivot

def clean_ontology_graph(uniprot_to_go, graph, ontologies):
    """
    Clean ontology-specific graphs by processing Uniprot mappings, layers, and nodes.

    Args:
        uniprot_to_go (dict): Mapping of UniProt IDs to GO terms.
        graph (networkx.DiGraph): The ontology graph.
        ontologies (list): Root nodes for ontologies.

    Returns:
        tuple: Cleaned graphs, updated leaves, selected nodes, max distances, and subgraphs.
    """
    try:
        selected_nodes, graph = map_uniprot_to_graph(uniprot_to_go, graph)

        graph_reverse = graph.reverse()
        max_dst, leaves, subgraphs, graph_reverse = get_leaves_by_ontology(ontologies, graph_reverse)

        graph = graph_reverse.reverse()
        graph = mark_lowest_nodes(ontologies, graph)
        clean_graphs, leaves, max_dst = get_clean_graph(leaves, max_dst, subgraphs, ontologies, graph)

        clean_graphs, leaves = clean_all_layers(clean_graphs, leaves, max_dst, ontologies)

        return clean_graphs, leaves, selected_nodes, max_dst, subgraphs
    except GraphCleaningDFSException as e:
        graph, ontologies = connect_pivot(ontologies, graph)
        return clean_ontology_graph(uniprot_to_go, graph, ontologies)

def process_cleaning(uniprot_to_go, graph, pivots, verbose=False):
    """
    Perform ontology graph cleaning.

    Args:
        uniprot_to_go (dict): Mapping of UniProt IDs to GO terms.
        graph (networkx.DiGraph): The ontology graph.
        pivots (list): Root nodes for ontologies.
        verbose (bool): print progress information.

    Returns:
        tuple: Cleaned graphs, updated pivots, max distances, and selected nodes.
    """
    try:
        clean_Graph, leaves, selected_nodes, max_dst, subgraph_ontology = clean_ontology_graph(uniprot_to_go, graph, pivots)
    except GraphCleaningDFSException as e:
        graph, pivots = connect_pivot(pivots, graph)
        clean_Graph, leaves, selected_nodes, max_dst, subgraph_ontology = clean_ontology_graph(uniprot_to_go, graph, pivots)

    if verbose == True:
        return clean_Graph,  pivots, leaves, selected_nodes, max_dst, subgraph_ontology
    else:
        return clean_Graph, pivots, max_dst, selected_nodes


def print_information_from_the_script(clean_Graph, pivot, verbose=False, selected_elements = None, most_distance = None, subgraph_ontology = None):
    for ontology in pivot:
        if len(clean_Graph[ontology].nodes) > 0 and 'name' in clean_Graph[ontology].nodes[ontology]:
            print("ontology: ", ontology, " - ", clean_Graph[ontology].nodes[ontology]['name'])
            print("ontology Graph information: ", clean_Graph[ontology])
            print("nodes in ontology "+ontology+": ",clean_Graph[ontology].nodes())
            print("edges in ontology "+ontology+": ",clean_Graph[ontology].edges())

        else:
            print("ontology: ", ontology)
            print("ontology Graph information: ", clean_Graph[ontology])
            print("node in ontology "+ontology+": ",clean_Graph[ontology].nodes())
            print("edges in ontology "+ontology+": ",clean_Graph[ontology].edges())


        if verbose == True:
            print('\n\n')
            print("selected element from add_UniprotID_to_node: ", selected_elements)
            print("most distance for all ontology: ", most_distance)
            print("subgraph for all ontology from add_layers_to_node: ", subgraph_ontology)

if __name__ == "__main__":
    Graph = obonet.read_obo("/home/pospim/Desktop/work/GOLizard/golizard/golizz/bash_app/meta/goslim_plant.obo")

    UniprotIDwithGOTerms = {'P05067': ['GO:0106003', 'GO:0045177', 'GO:0030424', 'GO:0009986', 'GO:0005911', 'GO:0035253', 'GO:0005905', 'GO:0030134', 'GO:0005737', 'GO:0005829', 'GO:0043198', 'GO:0043197', 'GO:0005769', 'GO:0005783', 'GO:0005788', 'GO:0005768', 'GO:0031904', 'GO:0070381', 'GO:0070062', 'GO:0005576', 'GO:0005615', 'GO:0005794', 'GO:0005796', 'GO:0005798', 'GO:0030426', 'GO:0034364', 'GO:0034363', 'GO:1990777', 'GO:0034362', 'GO:0016020', 'GO:0045121', 'GO:0005739', 'GO:0031594', 'GO:0005641', 'GO:0005634', 'GO:0043204', 'GO:0048471', 'GO:0005886', 'GO:0031093', 'GO:0048786', 'GO:0032991', 'GO:0043235', 'GO:0055037', 'GO:0005790', 'GO:0051233', 'GO:0045202', 'GO:0008021', 'GO:0032588', 'GO:0034361', 'GO:0030549', 'GO:0033130', 'GO:0097645', 'GO:0034185', 'GO:0042056', 'GO:0003682', 'GO:0003677', 'GO:0019899', 'GO:0046875', 'GO:0005109', 'GO:0001664', 'GO:1904399', 'GO:0043395', 'GO:0008201', 'GO:0042802', 'GO:0005158', 'GO:0005178', 'GO:0050750', 'GO:0046982', 'GO:0042803', 'GO:0051087', 'GO:0051425', 'GO:0050786', 'GO:0048018', 'GO:0000978', 'GO:0004867', 'GO:0030546', 'GO:0005102', 'GO:0046914', 'GO:0007189', 'GO:0007193', 'GO:0008344', 'GO:1990000', 'GO:0019731', 'GO:0019732', 'GO:0061844', 'GO:0008306', 'GO:0048143', 'GO:0002265', 'GO:0008088', 'GO:0016199', 'GO:0007409', 'GO:0019722', 'GO:0007155', 'GO:1904646', 'GO:0007417', 'GO:0008203', 'GO:0050890', 'GO:0048669', 'GO:0180011', 'GO:0050829', 'GO:0050830', 'GO:0016358', 'GO:0006897', 'GO:0030198', 'GO:0030900', 'GO:0007186', 'GO:0000086', 'GO:0045087', 'GO:0006878', 'GO:0035235', 'GO:0007612', 'GO:0007611', 'GO:0042157', 'GO:0007626', 'GO:0007617', 'GO:0007613', 'GO:0014005', 'GO:0001774', 'GO:0098815', 'GO:0006378', 'GO:1903523', 'GO:0090090', 'GO:0008285', 'GO:1902951', 'GO:0010629', 'GO:1900272', 'GO:1902894', 'GO:0010823', 'GO:0045665', 'GO:0010466', 'GO:1900181', 'GO:0000122', 'GO:0030178', 'GO:0050885', 'GO:0051402', 'GO:0070050', 'GO:0031175', 'GO:1990535', 'GO:0016322', 'GO:0007219', 'GO:0061903', 'GO:1905908', 'GO:0043065', 'GO:1902961', 'GO:0050867', 'GO:0032722', 'GO:0043280', 'GO:0007204', 'GO:0051091', 'GO:1904472', 'GO:0070374', 'GO:2000463', 'GO:2001238', 'GO:1904022', 'GO:0045745', 'GO:0010971', 'GO:0010628', 'GO:0045821', 'GO:0035066', 'GO:0050729', 'GO:0032731', 'GO:0032755', 'GO:0046330', 'GO:1900454', 'GO:1900273', 'GO:0043406', 'GO:0043410', 'GO:0051044', 'GO:0045931', 'GO:0090026', 'GO:0043525', 'GO:0045666', 'GO:0051092', 'GO:1901224', 'GO:0045429', 'GO:0033138', 'GO:0010800', 'GO:0042327', 'GO:0032092', 'GO:1904591', 'GO:0010739', 'GO:0051897', 'GO:0051247', 'GO:0001934', 'GO:0061098', 'GO:1900122', 'GO:1905898', 'GO:0032930', 'GO:2000406', 'GO:1902949', 'GO:0045944', 'GO:0032760', 'GO:0032729', 'GO:0051260', 'GO:0006468', 'GO:0051262', 'GO:0070206', 'GO:1903048', 'GO:1905906', 'GO:1900221', 'GO:1902950', 'GO:1903381', 'GO:0007176', 'GO:0010468', 'GO:0048169', 'GO:0043408', 'GO:0040014', 'GO:2000310', 'GO:0050730', 'GO:1905606', 'GO:0061097', 'GO:1905945', 'GO:0010469', 'GO:0150003', 'GO:0050803', 'GO:0034121', 'GO:0006357', 'GO:0006417', 'GO:0030111', 'GO:0070555', 'GO:0006979', 'GO:0001878', 'GO:0051563', 'GO:0001967', 'GO:0050808', 'GO:0051124', 'GO:0008542'],
             'P06213': ['GO:0030424', 'GO:0005901', 'GO:0032590', 'GO:0010008', 'GO:0009897', 'GO:0070062', 'GO:0005899', 'GO:0005770', 'GO:0005764', 'GO:0016020', 'GO:0032809', 'GO:0005886', 'GO:0043235', 'GO:0001540', 'GO:0005524', 'GO:0038024', 'GO:0005525', 'GO:0042802', 'GO:0043559', 'GO:0005009', 'GO:0043560', 'GO:0031994', 'GO:0031995', 'GO:0005159', 'GO:0043548', 'GO:0019904', 'GO:0004713', 'GO:0044877', 'GO:0051425', 'GO:0005198', 'GO:0032147', 'GO:0032148', 'GO:0030325', 'GO:0097242', 'GO:0005975', 'GO:0071363', 'GO:0032869', 'GO:0097062', 'GO:0008544', 'GO:0031017', 'GO:0007186', 'GO:0042593', 'GO:0003007', 'GO:0008286', 'GO:0007612', 'GO:0008584', 'GO:0030238', 'GO:0007613', 'GO:1990535', 'GO:0038083', 'GO:0018108', 'GO:0030335', 'GO:0008284', 'GO:0048639', 'GO:0045893', 'GO:0046326', 'GO:0045725', 'GO:0045821', 'GO:0033674', 'GO:0043406', 'GO:0043410', 'GO:0051446', 'GO:0045840', 'GO:0045429', 'GO:0014068', 'GO:0051897', 'GO:0001934', 'GO:0043243', 'GO:0002092', 'GO:0060267', 'GO:0046777', 'GO:0006468', 'GO:0031623', 'GO:0006898', 'GO:0006355', 'GO:0045995', 'GO:2000194', 'GO:0007169', 'GO:0150104', 'GO:0046718'],
             'A1A4S6': ['GO:0005829', 'GO:0048471', 'GO:0005886', 'GO:0005096', 'GO:0007010', 'GO:0043066', 'GO:0051056', 'GO:0007165']}

    pivot = ['GO:0003674', #molecular function
             'GO:0008150', #biological process
             'GO:0005575'] #cellular component

    #Verbose = False
    clean_Graph, pivot, most_distance, selected_elements = process_cleaning(UniprotIDwithGOTerms, Graph, pivot)
    print_information_from_the_script(clean_Graph, pivot)


    #Verbose = True
    #clean_Graph, pivot, leaves, selected_elements, most_distance, subgraph_ontology = returning_relevant_values(UniprotIDwithGOTerms, Graph, pivot, verbose=True)
    #print_information_from_the_script(clean_Graph, pivot, leaves = leaves, verbose=True, selected_elements = selected_elements, most_distance = most_distance, subgraph_ontology=subgraph_ontology)
