"""Produces graphs, tables, and saves node raw data to disk."""

import pandas as pd
from type_definitions import NodeType, OutputNodes
from pathlib import Path


def gen_dataframe(edit_tree, which_nodes: OutputNodes) -> pd.DataFrame:
    """Produces a DataFrame for the edit tree of a given guide node. This can then be saved to disk to minimise
    memory usage."""

    if which_nodes is OutputNodes.ALL:
        selected_nodes = [node for node in edit_tree.edit_nodes_all]
    elif which_nodes is OutputNodes.COMPLETE:
        selected_nodes = [node for node in edit_tree.edit_nodes_all if node.node_type is NodeType.COMPLETE]
    else:
        selected_nodes = [node for node in edit_tree.edit_nodes_all if node.progressed]

    output_values = [vars(node) for node in selected_nodes] # a list of dicts

    output_df = pd.DataFrame(output_values)
    output_df['sequence.seq'] = [node.sequence.seq for node in selected_nodes]

    return output_df


def save_edit_tree(guide_node, which_nodes: OutputNodes) -> Path:
    """Saves edit nodes and guide nodes to separate csv files. Can specify the output to include all edit nodes or
    only the complete, or progressed nodes."""

    df_of_edit_nodes = gen_dataframe(edit_tree=guide_node.edit_tree, which_nodes=which_nodes)

    tree = guide_node.guide_tree
    level = f'Level_{guide_node.guide_level:03d}'

    dest_addr = tree.log_path.parent / tree.id / level / guide_node.id / 'edit_node_values'
    dest_addr.parent.mkdir(parents=True, exist_ok=True)

    df_of_edit_nodes.to_csv(dest_addr)

    return dest_addr


def save_guide_tree(guide_tree) -> tuple:
    """Saves all guide node details to a DataFrame and outputs it as csv."""

    dest_addr = guide_tree.log_path.parent / Path(guide_tree.id) / Path('guide_node_values')
    dest_addr.parent.mkdir(parents=True, exist_ok=True)

    all_values = []
    for node in guide_tree.guide_nodes_all:
        vars_dict = vars(node)
        vars_dict['progressed_sequences'] = [tup[0].seq for tup in node.progressed_sequences]
        vars_dict['progressed_indices'] = [tup[1] for tup in node.progressed_sequences]
        vars_dict['guide']: node.guide_name.split('_')[1]
        if node.prior:
            vars_dict['used_priors'] = True
        else:
            vars_dict['used_priors'] = False
        if node.parent:
            vars_dict['parent_id'] = node.parent.id
        else:
            vars_dict['parent_id'] = None

        if vars_dict['progressed_indices']:
            vars_dict['end_mIndex'] = min(vars_dict['progressed_indices'])
        else:
            vars_dict['end_mIndex'] = 1

        all_values.append(vars_dict)

    guides_df = pd.DataFrame(all_values)

    guides_df.to_csv(dest_addr)

    return guides_df, dest_addr
