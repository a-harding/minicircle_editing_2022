"""Defines the class and associated functions for the guide node. This node represents a single guide RNA-mRNA
interaction and all the edits that this produces."""


import numpy as np
import uuid

import pandas as pd
import time

from docking import select_guides
from edit_trees import EditTree
from run_settings import editing_window, short_sequence_editing, max_no_grnas_subsequent, min_no_grnas_subsequent
from sequence_import import Sequence
from outputs_gen import save_edit_tree
from type_definitions import OutputNodes


class GuideNode:
    """Contains details of a single guide RNA interaction with an mRNA sequence."""

    def __init__(self, init_duplex: list, parent=None, seq_num=0, guide_num=0, guide_tree=None):
        self.parent = parent
        self.parents = [self.parent]
        self.children = []
        self.child_generation_investigated = False
        self.node_number: int

        if self.parent:
            self.guide_tree = self.parent.guide_tree
            self.prev_guides = self.parent.prev_guides + [self.parent.guide_name]
            self.guide_level = self.parent.guide_level + 1
        else:
            self.guide_tree = guide_tree
            self.prev_guides = []
            self.guide_level = 1
            self.node_number = 1

        self.seq_num = seq_num
        self.guide_num = guide_num
        self.id = f'{self.guide_tree.id}_L{self.guide_level}_S{seq_num}_G{guide_num}_{uuid.uuid4()}'

        self.duplex = init_duplex
        # initial_duplex is list with order [guide_name, mRNA_dock_index, mRNA_sequence, guide_index]
        self.guide_name = self.duplex[0]
        self.guide_sequence = self.guide_tree.guides_dict[self.guide_name]
        self.init_gIndex = self.duplex[3]

        self.init_sequence = self.duplex[2]
        self.init_dock_idx = self.duplex[1]
        self.init_mIndex = self.init_dock_idx + self.init_gIndex

        if short_sequence_editing:
            self.sequence_fragments = self.trim_editing_sequence()
            self.working_sequence = self.sequence_fragments[0]
            self.working_dock_idx = self.sequence_fragments[1]
            self.working_mIndex = self.sequence_fragments[2]
            self.leading_string = self.sequence_fragments[3]
            self.trailing_string = self.sequence_fragments[4]
        else:
            self.working_sequence = self.init_sequence
            self.working_dock_idx = self.init_dock_idx
            self.working_mIndex = self.init_mIndex

        # checks if the current duplex has been investigated prior and fetches the prior results if so
        self.prior = self.check_cache()
        # assigns the progressed sequences from the prior computation if applicable
        if self.prior:
            # print('*****\nPrior guide node computation used. No edit tree generated.\n*****')
            self.progressed_sequences = self.prior
        # in the absence of a prior computation, produces the edit tree as standard
        else:
            self.edit_tree = EditTree(guide_node=self)
            self.progressed_sequences = [(node.sequence, node.mIndex) for node in self.edit_tree.progressed_nodes]

        if self.progressed_sequences:
            self.is_terminal = False
        else:
            self.is_terminal = True

        if not self.prior:
            self.edit_tree_path = save_edit_tree(self, which_nodes=OutputNodes.ALL)
            self.cache_outputs()

        self.reassemble_sequence()
        self.errors_accumulated = self.total_errors()

    def gen_children(self):
        """Generates child guide nodes for each of the sequences that meet the criteria for progression."""

        for si, (sequence, mIndex) in enumerate(self.progressed_sequences):
            existing_children = self.guide_tree.get_existing_children(sequence.seq)
            if existing_children:
                # print('***\nChildren already exist\n***')
                self.children += existing_children
                # print('\n\nSleeping\n\n')
                # time.sleep(5)
                for child in existing_children:
                    child.parents.append(self)
                self.guide_tree.existing_children_uses += 1
            else:
                for gi, duplex in enumerate(select_guides(
                        messenger=sequence, previous_guides=self.prev_guides + [self.guide_name],
                        guides_dict=self.guide_tree.guides_dict, current_mIndex=mIndex - 5)):

                    # generate child nodes for the minimum number of guides to consider.
                    progressed_count = sum([1 for child in self.children if child.progressed_sequences])
                    conditions = (gi < min_no_grnas_subsequent) | ((gi < max_no_grnas_subsequent) & (progressed_count == 0))
                    if conditions:
                        new_guide_node = GuideNode(init_duplex=duplex, parent=self, seq_num=si + 1, guide_num=gi + 1)
                        self.children.append(new_guide_node)
                        new_guide_init_seq = pd.Series(data=new_guide_node.init_sequence.seq)
                        self.guide_tree.guide_nodes_series = self.guide_tree.guide_nodes_series.append(
                            new_guide_init_seq, ignore_index=True)
                        self.guide_tree.guide_nodes_all.append(new_guide_node)

        self.child_generation_investigated = True


    def cache_outputs(self) -> None:
        """Takes the key components of the guide node and places them into a results cache in the 'guide_tree'.
        The key components of the 'init_duplex' (containing initiating mRNA sequence, guide RNA, and associated
        indices), and output sequences (and indices) are stored in a pandas DataFrame.

        In the event that a previous guide_node has been created with the same initial mRNA sequence and guide RNA,
        then the output sequences and indices can be used from this prior determination, saving processing time."""

        var_dict = {
            'directory': str(self.edit_tree_path.parent),
            'guide_name': self.guide_name,
            'init_dock_index': self.duplex[1],
            'init_sequence': self.duplex[2].seq,
            'init_gIndex': self.duplex[3],
            'outputs': self.progressed_sequences
        }

        try:
            self.guide_tree.guide_node_cache = self.guide_tree.guide_node_cache.append(var_dict, ignore_index=True)
        except AttributeError:
            self.guide_tree.guide_node_cache = pd.DataFrame(var_dict)

        return None

    def check_cache(self) -> list:
        """Checks the cache of previously computed guide_nodes to see if the current duplex has been investigated
        already. If so, it returns a list containing the results from this prior node. If not, it returns
        an empty list."""

        try:
            df = self.guide_tree.guide_node_cache
        except AttributeError:
            return []

        cond_guide = df['guide_name'] == self.guide_name
        cond_dock = df['init_dock_index'] == self.duplex[1]
        cond_seq = df['init_sequence'] == self.duplex[2].seq
        cond_gindex = df['init_gIndex'] == self.duplex[3]

        condition_mask = cond_guide & cond_dock & cond_seq & cond_gindex

        matches = df[condition_mask]

        if len(matches):
            self.guide_tree.cache_uses += 1
            return matches['outputs'].iloc[0]
        else:
            return []

    def trim_editing_sequence(self) -> tuple:
        """Trims the mRNA sequence to only contain the sequence necessary for performing editing with this guide.
        It benefits from smaller memory size, processing and search times, and enables greater use of the guide_node
        cache, to remove unnecessary repetition."""

        init_seq_string = self.init_sequence.seq

        # the two indices that will be used to select out the mRNA region to be used in editing

        # start_idx uses the first base the guide docks at (init_dock_idx), with an additional n bases 3-prime
        # (as specified by 'editing_window') to allow for the 'CofoldType.EDITING_WINDOW' mode of RNA cofolding
        # if the docking site is too 3-prime to allow for the additional n bases, all 3-prime bases are included
        if (self.init_dock_idx - editing_window) < 0:
            start_idx = 0
            short_dock_idx = self.init_dock_idx
            short_mIndex = self.init_mIndex
        else:
            start_idx = self.init_dock_idx - editing_window
            short_dock_idx = editing_window
            short_mIndex = self.init_gIndex + editing_window

        # end_idx caps the length of the shortened mRNA. The selection of 50 additional bases is arbitrary but will
        # allow for the most extreme cases of consecutive deletions in a guide RNA's editing steps. The additional n
        # bases of the editing_window are included for the same reason as with start_idx
        end_idx = min(start_idx + self.guide_sequence.length + 50 + editing_window, self.init_sequence.length)

        leader_seq_string = init_seq_string[:start_idx]
        trailing_seq_string = init_seq_string[end_idx:]

        short_sequence_string = init_seq_string[start_idx:end_idx]
        short_sequence = Sequence(name=f'{self.init_sequence.name}_short',
                                  sequence=short_sequence_string,
                                  seq_is_5to3=False)

        return short_sequence, short_dock_idx, short_mIndex, leader_seq_string, trailing_seq_string

    def reassemble_sequence(self) -> None:
        """Re-attaches the leading and trailing mRNA sequences to the edited mRNA products if 'short_sequence_editing'
        is enabled. Also adjusts the mIndex to account for the additional 3-prime bases."""

        if not short_sequence_editing:
            return

        complete_progressed = []
        for si, (sequence, mIndex) in enumerate(self.progressed_sequences):
            complete_seq_string = f'{self.leading_string}{sequence.seq}{self.trailing_string}'
            complete_sequence = Sequence(name=sequence.name,
                                         sequence=complete_seq_string,
                                         seq_is_5to3=False)
            true_mIndex = mIndex + len(self.leading_string)
            complete_progressed.append((complete_sequence, true_mIndex))

        self.progressed_sequences = complete_progressed

        return

    def total_errors(self) -> list:
        """Counts the accumulated errors, from the 3-prime end of the mRNA up to the current mIndex, base-by-base,
        compared to the known edited mRNA sequence. Only used for QC of the model's development process."""

        known_edited = self.guide_tree.known_edited_sequence.seq
        total_errors = []
        for i, (sequence, index) in enumerate(self.progressed_sequences):
            current_errors = 0
            for known, model in zip(known_edited[:index], sequence.seq[:index]):
                if known != model:
                    current_errors += 1
            total_errors.append(current_errors)

        return total_errors
