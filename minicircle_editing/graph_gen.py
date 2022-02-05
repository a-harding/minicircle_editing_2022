"""Generates the graphviz outputs for guide trees and edit trees."""


import graphviz as gv
from pathlib import Path
from type_definitions import NodeType
from run_settings import maximum_edit_graph_nodes


def produce_match_string(edit_node) -> str:
    """Creates a combined string of pre- and post-edit sequences, paired with the guide sequence."""

    mismatch_set = {'gg', 'cc', 'aa', 'uu', 'ga', 'ag', 'ac', 'ca', 'cu', 'uc'}

    current_m_idx = edit_node.mIndex
    current_g_idx = edit_node.gIndex
    prior_m_seq = edit_node.init_sequence.seq
    post_m_seq = edit_node.sequence.seq
    g_seq = edit_node.guide.seq

    length_lead = 3
    length_trail = 4

    trimmed_prior_m_seq = f'{prior_m_seq[current_m_idx - length_lead:current_m_idx]}' \
                          f'{prior_m_seq[current_m_idx].upper()}' \
                          f'{prior_m_seq[current_m_idx + 1:current_m_idx + length_trail]}'
    trimmed_post_m_seq = f'{post_m_seq[current_m_idx - length_lead:current_m_idx]}' \
                         f'{post_m_seq[current_m_idx].upper()}' \
                         f'{post_m_seq[current_m_idx + 1:current_m_idx + length_trail]}'
    trimmed_g_seq = f'{g_seq[current_g_idx - length_lead:current_g_idx]}' \
                    f'{g_seq[current_g_idx].upper()}' \
                    f'{g_seq[current_g_idx + 1:current_g_idx + length_trail]}'

    if len(trimmed_g_seq) < len(trimmed_post_m_seq):
        missing = len(trimmed_post_m_seq) - len(trimmed_g_seq)
        filler = missing * ' '
        trimmed_g_seq = f'{trimmed_g_seq}{filler}'

    def compare_sequences(s1, s2):
        match_string = ''
        for b1, b2 in zip(s1, s2):
            if b2 == ' ':
                result = ' '
            else:
                pair = (b1 + b2).lower()
                if pair in mismatch_set:
                    result = '.'
                else:
                    if pair in {'gu', 'ug'}:
                        result = 'o'
                    else:
                        result = '|'
            match_string = match_string + result
        return match_string

    prior_match = compare_sequences(trimmed_prior_m_seq, trimmed_g_seq)
    post_match = compare_sequences(trimmed_post_m_seq, trimmed_g_seq)

    compound_string = f'{trimmed_prior_m_seq}\n{prior_match}\n{trimmed_g_seq}\n{post_match}\n{trimmed_post_m_seq}'

    return compound_string


def graph_edit_tree(edit_tree) -> gv.Digraph:
    """Produces dot format Diagraph from the input list of edit nodes."""

    guide_tree = edit_tree.guide_node.guide_tree
    level = f'Level_{edit_tree.guide_node.guide_level:03d}'

    file_name = Path(f'edit_tree_{edit_tree.id}')
    dest_addr = guide_tree.log_path.parent / guide_tree.id / level / edit_tree.guide_node.id / file_name
    dest_addr.parent.mkdir(parents=True, exist_ok=True)
    graph = gv.Digraph(engine='dot', format='svg')

    all_nodes = edit_tree.edit_nodes_all

    if len(all_nodes) > maximum_edit_graph_nodes:
        print('Edit node total exceeds acceptable graphable threshold. No edit tree graph generated.')
        return None

    def create_edit_node_text(edit_node) -> str:
        """Produces text content for nodes in the tree graph."""

        # if edit_node.node_type is NodeType.ROOT:
        dock = f'Dock: {edit_node.edit_tree.guide_node.init_dock_idx}'
        id = f'ID: {edit_node.id}'
        action_log = f'Action log: {edit_node.action_log[-5:]}'
        type = f'Node Type: {edit_node.node_type}'
        gIdx = f'Guide Index: {edit_node.gIndex}'
        mIdx = f'Seq Index: {edit_node.mIndex}'
        action = f'Action: {edit_node.action}'
        try:
            mfe = f'MFE: {edit_node.mfe: 0.5f}'
            align = edit_node.alignment
            prob = f'Prob {edit_node.probability:0.3f}'
            prob_prod = f'Prob_prod: {edit_node.probability_product:0.3f}'
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
        match = produce_match_string(edit_node=edit_node)

        if edit_node.parent:
            text = f'{mismatches}\n{prob}\n{prob_prod}\n{match}'
        else:
            text = 'Root Node'

        return text

    colourschemes = {0: 'bupu', 1: 'purd', 2: 'pubu', 3: 'bugn', 4: 'blues', 5: 'greys', 6: 'oranges',
                     7: 'purples', 8: 'reds'}
    scheme = colourschemes[0]
    graph.attr('node', style='filled', fontname='monospace', colorscheme=scheme + str(2 + 2))

    for node in all_nodes:
        node_text = create_edit_node_text(node)
        # sets colour of node by its number of mismatches
        if node.progressed:
            node_shape = 'square'
        else:
            # node_shape = 'ellipse'
            node_shape = 'circle'

        # colour_num = min(1, node.probability_product) * 100
        # colour = f'grey{colour_num:0.0f}'
        # if colour_num < 60:
        #     text_colour = 'white'
        # else:
        #     text_colour = 'black'
        # graph.attr('node', shape=node_shape, fontcolor = text_colour, fillcolor = colour)
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
            actions = {'p': 'blue', 'd': 'red', 'i': 'green'}
            colour = actions[child.action]
            graph.edge(tail_name=parent.id, head_name=child.id, color=colour)

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

        id_split = guide_node.id.split('-')[0].split('_')
        id = f'{id_split[0]}\n{id_split[1]}\n{id_split[4]}'
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

        text = f'{id}\n{dock_idx}\n{tot_err}'

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

        if len(node.parents) > 1:
            colour = 'grey94'
        else:
            colour = 'grey99'

        # if len(node.errors_accumulated):
        #     if min(node.errors_accumulated) == 0:
        #         colour = 'green'
        # else:
        #     colour = 'red'

        if node.correct_dock:
            colour = 'green'

        graph.attr('node', shape='ellipse', fillcolor=colour)
        graph.node(str(node.id), node_text)

    # adds edge for each parent-child pair, even if the child has multiple parents
    for node in all_nodes:
        # if node.parent:
        shapes = {1: 'normal', 2: 'box', 3: 'curve', 4: 'dot', 5: 'inv', 6: 'tee'}
        colors = {1: 'green', 2: 'gold', 3: 'blue', 4: 'hotpink', 5: 'black', 6: 'azure4'}

        seq_num = node.seq_num

        for child in node.children:
            guide_num = child.guide_num
            graph.attr('edge', arrowhead=shapes[1], color=colors[min(guide_num, 6)])
            # tail_name = parent, head_name = child
            graph.edge(label=child.g_name_short,
                       tail_name=node.id,
                       head_name=child.id)

    graph.attr('edge', style='solid')

    graph.render(str(dest_addr), view=False)

    return graph
