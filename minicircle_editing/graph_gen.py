"""Generates the graphviz outputs for guide trees and edit trees."""


import graphviz as gv
from pathlib import Path
from type_definitions import NodeType
from run_settings import maximum_edit_graph_nodes


def graph_edit_tree(edit_tree) -> gv.Digraph:
    """Produces dot format Diagraph from the input list of edit nodes."""

    file_name = Path(f'edit_tree_{edit_tree.id}')
    dest_addr = edit_tree.guide_node.guide_tree.log_path.parent / Path(edit_tree.guide_node.id) / file_name
    dest_addr.parent.mkdir(parents=True, exist_ok=True)
    graph = gv.Digraph(engine='dot', format='svg')

    all_nodes = edit_tree.edit_nodes_all

    if len(all_nodes) > maximum_edit_graph_nodes:
        print('Edit node total exceeds acceptable graphable threshold. No edit tree graph generated.')
        return None

    colourschemes = {0: 'bupu', 1: 'purd', 2: 'pubu', 3: 'bugn', 4: 'blues', 5: 'greys', 6: 'oranges',
                     7: 'purples', 8: 'reds'}
    scheme = colourschemes[0]
    graph.attr('node', style='filled', colorscheme=scheme + str(2 + 2))

    def create_edit_node_text(edit_node) -> str:
        """Produces text content for nodes in the tree graph."""

        id = f'ID: {edit_node.id}'
        action_log = f'Action log: {edit_node.action_log}'
        type = f'Node Type: {edit_node.node_type}'
        gIdx = f'Guide Index: {edit_node.gIndex}'
        mIdx = f'Seq Index: {edit_node.mIndex}'
        action = f'Action: {edit_node.action}'
        try:
            mfe = f'MFE: {edit_node.mfe: 0.5f}'
            align = edit_node.alignment
            prob = f'Probability {edit_node.probability:0.4f}'
            prob_prod = f'Probability (Product): {edit_node.probability_product:0.4f}'
        except AttributeError:
            mfe = f'MFE not calculated'
            align = 'No alignment'
            prob = f'Probability not calculated'
            prob_prod = f'Probability (Product) not calculated'
        mismatches = f'Mismatches: {edit_node.mismatches}'
        # mismatches_full = f'Mismatches (f): {node.mismatches_full}'
        seq_pre = f'Seq init: {edit_node.init_sequence.seq[edit_node.mIndex - 5:edit_node.mIndex + 5]}'
        seq_post = f'Seq pos: {edit_node.sequence.seq[edit_node.mIndex - 5:edit_node.mIndex + 5]}'
        guide = f'gSeq: {edit_node.guide.seq[edit_node.edit_tree.init_gIndex:]}'
        # pairs = f'Pairs: {node.pairs}'

        text = f'{type}\n{mIdx} {gIdx}\n{action_log}\n{mfe}\n{mismatches}\n{seq_pre}\n{seq_post}\n{prob}\n{prob_prod}'

        return text

    for node in all_nodes:
        node_text = create_edit_node_text(node)
        # sets colour of node by its number of mismatches
        if node.progressed:
            node_shape = 'square'
        else:
            node_shape = 'ellipse'
        graph.attr('node', shape=node_shape, fillcolor=str(min(int(node.mismatches), 2 + 1) + 1))
        graph.node(str(node.id), node_text)

    # # adds edge for each parent-child pair, even if the child has multiple parents
    # for node in all_nodes:
    #     if node.node_type is NodeType.ROOT:
    #         continue
    #     else:
    #         for parent in node.parents:
    #             # tail_name = parent, head_name = child
    #             graph.edge(tail_name=parent.action_log, head_name=node.action_log)


    for parent in all_nodes:
        # tail_name = parent, head_name = child
        for child in parent.children:
            graph.edge(tail_name=parent.id, head_name=child.id)

    graph.attr('edge', style='solid')

    graph.render(str(dest_addr), view=False)

    return graph

def graph_guide_tree(guide_tree):
    """Produces dot format Diagraph from the input list of nodes. Outputs a txt file detailing the and the
     node ediges and a pdf containing the graph itself."""

    dest_addr = guide_tree.log_path.parent / Path(guide_tree.id) / Path('guide_tree')
    dest_addr.parent.mkdir(parents=True, exist_ok=True)
    graph = gv.Digraph(engine='dot')

    all_nodes = guide_tree.guide_nodes_all

    def create_guide_node_text(guide_node) -> str:
        """Produces text content for nodes in the guide tree graph."""

        init_dock_idx = guide_node.duplex[1]
        init_sequence = guide_node.duplex[2]
        prog_sequences_dict = guide_node.progressed_sequences
        if prog_sequences_dict:
            prog_count = len(prog_sequences_dict)
        else:
            prog_count = 'None'
        prog_seqs = ''
        if prog_sequences_dict:
            for i, (seq, idx) in enumerate(prog_sequences_dict):
                prog_seqs = prog_seqs + f'{i}: {seq.seq[:idx]}\n'
            else:
                prog_seqs = 'None'

        id = f'ID: {guide_node.id}'
        # type = f'Node Type: {guide_node.node_type}'
        # gIdx = f'Guide Index: {guide_node.gIndex}'
        dock_idx = f'Dock Index: {init_dock_idx}'
        # mIdx = f'Seq Index: {init_mIndex}'
        guide = f'{guide_node.guide_name}'
        progressed_count = f'Progressed: {prog_count}'
        progressed_seqs = f'Progressed sequences:\n{prog_seqs}'
        # action = f'Action: {guide_node.action}'
        # prob_prod = f'Probability: {node.probability_product:0.4f}'
        # mismatches = f'Mismatches: {guide_node.mismatches}'
        # mismatches_full = f'Mismatches (f): {node.mismatches_full}'
        seq_pre = f'Seq init: {guide_node.duplex[2].seq[init_dock_idx:init_dock_idx + 40]}'
        # seq_post = f'Seq pos: {guide_node.sequence.seq[15:30]}'
        # guide = f'gSeq: {guide_node.guide.seq[guide_node.edit_tree.init_gIndex:]}'
        tot_err = f'Errors in progressed: {guide_node.errors_accumulated}'

        text = f'{guide}\n{dock_idx}\n{progressed_count}\n{tot_err}'

        return text

    graph.attr('node', style='filled')
    for node in all_nodes:
        init_dock_idx = node.duplex[1]
        init_sequence = node.duplex[2]
        node_text = create_guide_node_text(node)
        # sets colour of node by progression across the mRNA. Uses init_mIndex compared with current mRNA length
        # to define a percentage of progression.
        progress = ((init_dock_idx) / init_sequence.length) * 100
        # graph.attr('node', shape='ellipse', fillcolor=f'gray{100 - progress: 0.0f}')
        graph.attr('node', shape='ellipse', fillcolor=f'gray99')
        graph.node(str(node.id[-20:]), node_text)

    # adds edge for each parent-child pair, even if the child has multiple parents
    for node in all_nodes:
        if node.parent:
            shapes = {1: 'normal', 2: 'box', 3: 'curve', 4: 'dot', 5: 'inv', 6: 'tee'}
            colors = {1: 'green', 2: 'gold', 3: 'blue', 4: 'hotpink', 5: 'black', 6: 'azure4'}
            guide_num = node.guide_num
            seq_num = node.seq_num
            graph.attr('edge', arrowhead=shapes[guide_num], color=colors[seq_num])
            # tail_name = parent, head_name = child
            graph.edge(tail_name=node.parent.id[-20:], head_name=node.id[-20:])

    graph.attr('edge', style='solid')

    graph.render(str(dest_addr), view=False)

    return graph
