#!/usr/bin/python3
import obonet
import networkx as nx
import time

EDGE_COLORS = {
    'is_a': '#000000',                    # black
    'part_of': '#000080',                 # dark blue
    'regulates': '#FFD700',               # gold
    'positively_regulates': '#00FF00',    # green
    'negatively_regulates': '#FF0000',    # red
    'occurs_in': '#40E0D0',               # turquoise
    'capable_of': '#ADD8E6',              # light blue
    'capable_of_part_of': '#FFA500',      # orange
    'has_part': '#FFC0CB',                # pink
    'default_color': '#FFFFFF'            # white
}
NODE_COLORS = {
    'GO:0003674': '#1f77b4',  # Blue for molecular function
    'GO:0008150': '#2ca02c',  # Green for biological process
    'GO:0005575': '#ff7f0e',  # Orange for cellular component
}
file_timestamp = int(time.time())

def color_nodes_by_namespace(agraph, graph):
    """
    Add colors to nodes based on their GO namespace.

    Args:
        agraph (pgv.AGraph): PyGraphviz AGraph object.
        graph (networkx.Graph): Original GO graph.

    Returns:
        pgv.AGraph: Updated AGraph with colored nodes.
    """
    for node in graph.nodes():
        node_agraph = agraph.get_node('GO_TERM'+node)
        node_agraph.attr['style'] = 'filled'
        node_agraph.attr["fontcolor"] = "white"
        namespace = graph.nodes[node].get('namespace', '')

        if namespace == 'molecular_function':
            color = NODE_COLORS['GO:0003674']
        elif namespace == 'biological_process':
            color = NODE_COLORS['GO:0008150']
        elif namespace == 'cellular_component':
            color = NODE_COLORS['GO:0005575']
        else:
            color = 'black'

        node_agraph.attr['color'] = color
        node_agraph.attr['fillcolor'] = color

    return agraph


def split_label(name, max_line_length=20):
    """
    Split a label into multiple lines for better readability.

    Args:
        name (str): Node name.
        max_line_length (int): Maximum characters per line.

    Returns:
        str: Multiline label.
    """
    words = name.replace('_', ' ').split()
    lines = []
    current_line = []

    for word in words:
        if sum(len(w) for w in current_line) + len(current_line) + len(word) <= max_line_length:
            current_line.append(word)
        else:
            lines.append(" ".join(current_line))
            current_line = [word]

    if current_line:
        lines.append(" ".join(current_line))

    return "\n".join(lines)

def get_node_label(node, graph):
    """
    Generate a label for a node based on its attributes.

    Args:
        node (str): Node ID.
        graph (networkx.Graph): GO graph.

    Returns:
        str: Node label.
    """
    if node.startswith('GO_TERM'):
        return node[7:]
    if node.startswith('name_GO_TERM'):
        return split_label(graph.nodes[node[12:]].get('name', ''))
    if node.startswith('UIDs'):
        uniprot_ids = graph.nodes[node[4:]].get('UniprotID', [])
        if isinstance(uniprot_ids, str):
            return uniprot_ids
        if len(uniprot_ids) <= 6:
            return ",\n".join(uniprot_ids)
        return f"Number of proteins referencing this GO Term: {len(uniprot_ids)}"
    return "Unknown"

def label_edges(graph):
    """
    Assign edge types to a graph based on node relationships.

    Args:
        graph (networkx.Graph): GO graph.

    Returns:
        networkx.Graph: Graph with assigned edge types.
    """
    for src, dst, data in graph.edges(data=True):
        relationships = graph.nodes[src].get('relationship', [])
        for rel in relationships:
            rel_type, target = rel.split()
            if target == dst:
                data['type_of_edge'] = rel_type
                break
        else:
            if dst in graph.nodes[src].get('is_a', []):
                data['type_of_edge'] = 'is_a'
            else:
                data['type_of_edge'] = 'default_color'
    return graph


def color_edges(agraph, graph):

    for edge in graph.edges(data=True):
        e_agraph = None
        if 'UniprotID' in graph.nodes[edge[0]].keys() and graph.nodes[edge[0]]['UniprotID'] != "lowest":
            e_agraph = agraph.get_edge("UIDs"+edge[0], "GO_TERM"+edge[1])
        else:
            e_agraph = agraph.get_edge("name_GO_TERM"+edge[0], "GO_TERM"+edge[1])

        if "type_of_edge" in edge[2]:
            e_agraph.attr['color'] = EDGE_COLORS[edge[2]["type_of_edge"]]
        else:
            e_agraph.attr['color'] = EDGE_COLORS['default_color']

    return agraph


def generate_graphviz(graph):
    """
    Generate a Graphviz visualization of the GO graph.

    Args:
        Graph (networkx.Graph): GO graph.

    Returns:
        pgv.AGraph: Graphviz AGraph object.
    """
    label_edges(graph)
    if "new_conecting_node" in graph.nodes():
        graph.remove_node("new_conecting_node")

    graph_copy = graph.copy()

    for node in list(graph_copy.nodes()):
        subgraph = nx.DiGraph()
        subgraph.add_node('GO_TERM'+node)
        subgraph.add_node('name_GO_TERM'+node)

        if 'UniprotID' in graph.nodes[node] and graph.nodes[node]['UniprotID'] != "lowest":
            subgraph.add_node('UIDs' + node)
            subgraph.add_edge('name_GO_TERM' + node, 'UIDs' + node)

        subgraph.add_edge('GO_TERM' + node, 'name_GO_TERM' + node)
        preds = list(graph.predecessors(node))
        succs = list(graph.successors(node))
        graph.remove_node(node)

        for subnode in subgraph.nodes():
            graph.add_node(subnode)
        for subedge in subgraph.edges():
            graph.add_edge(*subedge)

        for pred in preds:
            graph.add_edge(pred, 'GO_TERM'+node)

        for succ in succs:
            target = 'UIDs' + node if 'UniprotID' in graph_copy.nodes[node] and graph_copy.nodes[node]['UniprotID'] != "lowest" else 'name_GO_TERM' + node
            graph.add_edge(target, succ)

    agraph = nx.nx_agraph.to_agraph(graph)
    agraph = color_edges(agraph, graph_copy)
    agraph = color_nodes_by_namespace(agraph, graph_copy)

    for node in list(graph_copy.nodes()):
        cluster_name = f'cluster_{node}'
        subgraph_nodes = ['GO_TERM'+node, 'name_GO_TERM'+node]

        if 'UniprotID' in graph_copy.nodes[node].keys() and graph_copy.nodes[node]['UniprotID'] != "lowest":
            subgraph_nodes.append('UIDs'+node)

        agraph.add_subgraph(subgraph_nodes, name=cluster_name, shape='box')

        for node in subgraph_nodes:
            agraph.get_node(node).attr['shape'] = 'box'
            agraph.get_node(node).attr['label'] = get_node_label(node, graph_copy)

        for i in range(len(subgraph_nodes) - 1):
            if agraph.has_edge(subgraph_nodes[i], subgraph_nodes[i + 1]):
                agraph.get_edge(subgraph_nodes[i], subgraph_nodes[i + 1]).attr['style'] = 'invis'

    agraph.graph_attr['rankdir'] = 'TB'

    return agraph



if __name__ == "__main__":
    go_graph_file = "./input.obo"
    Graph = obonet.read_obo(go_graph_file)

    for node in Graph.nodes():
        Graph.nodes[node]["Layer_ofGO_0003674"] = int(Graph.nodes[node]["Layer_ofGO_0003674"][0])
        if 'UniprotID' in Graph.nodes[node]:
            Graph.nodes[node]['UniprotID'].extend(["extra1", "extra2", "extra3", "extra4", "extra5"])

    agraph = generate_graphviz(Graph)
    output_file = f"./dotfiles/output_graphviz_ontology_{file_timestamp}.dot"
    print(f'Output file: {output_file}')
    agraph.write(output_file)