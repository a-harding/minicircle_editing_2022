"""Builds trees of the entire editing process for a single mRNA. Each node represents a single guide RNA binding
and all the edits of that guide."""


import pandas as pd
import time
import pickle
from log_gen import convert_timestamp, append_log
from guide_node import GuideNode
from pathlib import Path
from graph_gen import graph_guide_tree
from outputs_gen import save_guide_tree


class GuideTree:
    """Container for all GuideNode objects in a single initial guide binding to the unedited mRNA sequence."""

    def __init__(self, initial_duplex, guides_dict, gene, log_path, edited_seq):
        self.gene = gene
        self.id = f'{self.gene}_{initial_duplex[0]}_mDI{initial_duplex[1]}_gI{initial_duplex[3]}'
        # initial_duplex is list with order [guide_name, mRNA_dock_index, mRNA_sequence, guide_index]
        self.guides_dict = guides_dict

        self.known_edited_sequence = edited_seq

        self.log_path = log_path
        self.timestamp = time.perf_counter()
        self.log(message=[f'\nNew guide tree: {self.id}\n'])

        self.guide_node_cache: pd.DataFrame
        self.cache_uses = 0

        self.root = GuideNode(guide_tree=self, init_duplex=initial_duplex)

        self.guide_nodes_all = [self.root]
        self.guide_nodes_current = [self.root]

        self.guide_levels = 1

        self.is_complete = False

        self.log(prev_nodes=0)

        while (not self.is_complete) and (self.guide_levels < 30):
            nodes_to_process = len(self.guide_nodes_current)
            self.grow_tree()
            self.guide_levels += 1
            self.log(prev_nodes=nodes_to_process)
            if not self.guide_nodes_current:
                self.is_complete = True

        print(f'***********************************************************\n{self.cache_uses} cache uses.')
        self.graph = graph_guide_tree(self)
        self.output_data = save_guide_tree(self)

    def log(self, prev_nodes=0, message=None) -> None:
        """Appends the latest guide tree-building event to the log file."""

        if message:
            append_log(log_path=self.log_path, new_events=message)
        else:
            last_timestamp = self.timestamp
            self.timestamp = time.perf_counter()
            hrs, mins, secs = convert_timestamp(time_start=last_timestamp, time_end=self.timestamp)

            lvl_complete_msg = f'Level {self.guide_levels - 1} complete. Nodes processed: {prev_nodes}\n'
            elapsed = f'Elapsed time: {hrs} hours, {mins} minutes, {secs:0.2f} seconds\n'
            gen_msg = f'Level {self.guide_levels} guide nodes generated. Nodes created: {len(self.guide_nodes_current)}\n'

            append_log(log_path=self.log_path, new_events=[lvl_complete_msg, elapsed, gen_msg])

        return None

    def grow_tree(self):
        """Grows the tree by instantiating child nodes for all the non-terminal guide nodes."""

        next_nodes = []
        nodes_num_to_process = len(self.guide_nodes_current)
        for guide_node in self.guide_nodes_current:
            print(f'Processing parent guide node {guide_node.node_number} of {nodes_num_to_process}')
            if not guide_node.is_terminal:
                guide_node.gen_children()
                for i, child_node in enumerate(guide_node.children):
                    child_node.node_number = len(next_nodes) + 1 + i
                    print(f'Guide node Level {child_node.guide_level}, Node {child_node.node_number} created.')
                next_nodes += guide_node.children
        # self.save_guide_nodes()
        self.guide_nodes_all += next_nodes
        self.guide_nodes_current = next_nodes

        return None

    def save_guide_nodes(self):
        """Pickles guide nodes from previous levels, and saves them to a file. Freeing up memory and enabling
        loading at a later stage."""

        for guide_node in self.guide_nodes_all:
            if isinstance(guide_node, GuideNode):
                file_path = self.log_path.parent / Path(guide_node.id) / Path('guide_node')
                file_path.parent.mkdir(parents=True, exist_ok=True)
                with open(file_path, mode='wb') as f:
                    pickle.dump(guide_node, f)
                guide_node = file_path
