#!/usr/bin/python3

import re

def add_links_to_go(svg_content):
    """
    Add clickable links to GO terms in an SVG file.

    Args:
        svg_content (str): Content of the SVG file.

    Returns:
        str: Modified SVG content with links added to GO terms.
    """
    # Regular expression to match the GO term group
    pattern = re.compile(
        r'(<g id="node\d+" class="node">\n<title>GO_TERM(GO:\d+)</title>(.*?)</g>)',
        re.DOTALL
    )

    def replace_link(match):
        group_content = match.group(1)
        go_term = match.group(2)
        return f'<a xlink:href="https://www.ebi.ac.uk/QuickGO/term/{go_term}">{group_content}</a>'

    new_svg_content = re.sub(pattern, replace_link, svg_content)
    return new_svg_content

def process_svg_file(file_path):
    """
    Process an SVG file to add clickable links to GO terms.

    Args:
        file_path (str): Path to the SVG file.
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        svg_content = file.read()

    updated_svg = add_links_to_go(svg_content)

    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(updated_svg)

if __name__ == "__main__":

    svg_file_path = './output_graphviz_new_conecting_node_1720103224.svg'
    process_svg_file(svg_file_path)
