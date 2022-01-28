"""Produces trees for every edit of a single guide RNA. Each node represents an edit decision, whereby the editing
machinery either deletes or inserts a U, or proceeds to the subsequent base without editing."""
import time

import RNA
import numpy as np
import pandas as pd
from edit_node import EditNodeRoot, EditNodeChild
from run_settings import bulk_cofold, sequences_to_progress, min_mfe_to_progress
from type_definitions import NodeType
from graph_gen import graph_edit_tree
from sequence_import import Sequence
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import os
from functools import lru_cache


class EditTree:
    """Contains all edit nodes for a single guide RNA's series of edits."""

    def __init__(self, guide_node):
        # print('Edit tree initiated.')
        self.guide_node = guide_node

        self.id = f'{self.guide_node.id}_ET'

        self.guide_sequence = self.guide_node.guide_sequence
        self.init_gIndex = self.guide_node.init_gIndex

        self.init_sequence = self.guide_node.working_sequence
        self.init_dock_idx = self.guide_node.working_dock_idx
        self.init_mIndex = self.guide_node.working_mIndex

        self.root = EditNodeRoot(edit_tree=self, init_seq=self.init_sequence, action='r')
        self.edit_nodes_all = [self.root]
        self.edit_nodes_next = []

        self.edit_levels = 1

        self.is_complete = False

        # print('sleeping')
        # time.sleep(5)

        while not self.is_complete:
            self.grow_tree()

        if bulk_cofold:
            self.calc_probabilities(self.edit_nodes_all)

        # self.candidates = [edit_node for edit_node in self.edit_nodes_all if edit_node.node_type is NodeType.COMPLETE]

        self.progressed_nodes = self.select_progressed_nodes()
        self.graph = graph_edit_tree(self)

    def calc_probabilities(self, nodes_list) -> None:
        """Calls the cofold_sequence function for each node provided, then calculates the probabilities for all."""

        k = 1.986e-3  # boltzmann constant, kcal/(mol.K)
        T = 310  # mammalian body temperature, K

        if bulk_cofold:
            working_nodes = [edit_node for edit_node in nodes_list if edit_node.node_type is NodeType.COMPLETE]
        else:
            working_nodes = sorted(nodes_list)

        # print(f'Number of nodes to cofold: {len(working_nodes)}')

        if not working_nodes:
            return None

        # # multiprocessing option 1
        # with ProcessPoolExecutor() as executor:
        #     working_nodes = executor.map(cofold_sequence_func, working_nodes)
        # mfes = [edit_node.mfe for edit_node in working_nodes]

        if len(working_nodes) > 100:
            # multiprocessing option 2
            mfes = []
            cofold_sequences = [node.cofold_string for node in working_nodes]
            with mp.Pool(processes=os.cpu_count() - 1) as p:
                cofold_results = p.map(RNA.cofold, cofold_sequences)

                # mfe_min = min([result[1] for result in cofold_results])

                for node, cofold_result in zip(working_nodes, cofold_results):
                    node.alignment, node.mfe = cofold_result
                    mfes.append((node.mfe))
        else:
            # single core option
            mfes = [edit_node.cofold_sequence() for edit_node in working_nodes]

        mfe_min = min(mfes)

        @lru_cache()
        def max_downstream_prob(node) -> float:
            """Recursive determination of the highest downstream 'probability' result."""

            best_probability = node.probability
            for child in node.children:
                child_probability = max_downstream_prob(child)
                if child_probability > best_probability:
                    best_probability = child_probability

            return best_probability

        # Option 1
        # Calculate the probability for all nodes of the current level
        for edit_node in working_nodes:
            edit_node.probability = np.exp((mfe_min - edit_node.mfe) / (k * T))

        # split out the 'probability_product' calculation to enable the upgrading of 'probability' value
        # if a node has a downstream node with a higher 'probability' result
        for edit_node in working_nodes:
            edit_node.probability = max_downstream_prob(edit_node)
            if bulk_cofold:
                edit_node.probability_product = edit_node.probability
            else:
                edit_node.probability_product = edit_node.parent.probability_product * edit_node.probability
                edit_node.probability_calculated = True

                if edit_node.node_type is NodeType.MERGED:
                    number_merged = sum([sibling.probability_calculated for sibling in edit_node.merge_siblings])
                    if number_merged == len(edit_node.merge_siblings):
                        sum_prob = sum([sibling.probability_product for sibling in edit_node.merge_siblings])
                        for sibling in edit_node.merge_siblings:
                            sibling.probability_product = sum_prob

        # for edit_node in sorted(working_nodes):
        #     if edit_node.node_type is NodeType.MERGED:
        #         if not edit_node.merge_complete:
        #             sum_prob = sum([sibling.probability_product for sibling in edit_node.merge_siblings])
        #             for sibling in edit_node.merge_siblings:
        #                 sibling.probability_product = sum_prob
        #                 sibling.merge_complete = True

        for edit_node in working_nodes:
            edit_node.reset_node_type()

            # prob_prod_total = 0
            #         for parent in edit_node.parents:
            #             prob_prod_total += parent.probability_product
            #         edit_node.combined_parent_prob_prod = prob_prod_total
            #         edit_node.probability_product = edit_node.combined_parent_prob_prod * edit_node.probability
            #     else:
            #         edit_node.probability_product = edit_node.parent.probability_product * edit_node.probability
            #     edit_node.reset_node_type()


            # if edit_node.node_type is NodeType.ROOT:
            #     print('Sleeping...')
            #     time.sleep(5)
            #     print('Done sleeping.')

        # # Option 2
        # mfes np.array(mfes)
        # parent_prob_prods = np.array([edit_node.parent.probability_product for edit_node in self.edit_nodes_next])
        # node_probs = np.exp((mfe_min - mfes) / (k * T))
        # node_prob_prods = parent_prob_prods * node_probs
        # for i, edit_node in enumerate(self.edit_nodes_next):
        #     edit_node.probability = node_probs[i]
        #     edit_node.probability_product = node_prob_prods[i]
        #     edit_node.reset_node_type()

        return None

    # def merge_nodes(self, node1, node2) -> EditNodeChild:
    #     """Merge two edit nodes with the same sequence generated by different edit pathways."""
    #
    #     print('Nodes merged')
    #
    #     # do merge
    #
    #     # OPTION 1 - merge into one of the existing nodes
    #     if node1.merged or not node2.merged:
    #         merge_node = node1
    #         donor_node = node2
    #     else:
    #         merge_node = node2
    #         donor_node = node1
    #
    #     merge_node.parents += donor_node.parents
    #     merge_node.action = f'{merge_node.action}{donor_node.action}'
    #     merge_node.id = f'{merge_node.parent.id}{merge_node.action}'
    #     merge_node.action_log = f'{merge_node.parent.action_log}({merge_node.action})'
    #
    #     for parent in donor_node.parents:
    #         if donor_node.action[0] == 'd':
    #             parent.child_deledit[0] = merge_node
    #             parent.children[1] = merge_node
    #         elif donor_node.action[0] == 'i':
    #             parent.children_nondel[1] = merge_node
    #             parent.children[1] = merge_node
    #         else:
    #             parent.children_nondel[0] = merge_node
    #             parent.children[0] = merge_node
    #
    #     merge_node.merged = True
    #
    #     # OPTION 2 - make a new node of the same level that treats the donor nodes as parents
    #
    #     # edit_level = 0
    #     # actions = ''
    #     # for node in node_list:
    #     #     actions = f'{actions}{node.action}'
    #     #     if node.edit_level < edit_level:
    #     #         edit_level = node.edit_level
    #
    #     return merge_node

    def get_match_indices(self, node_list, search_item):
        """Gets a list of indices where the item in a list matches an item being searched for."""

        indices_list = []
        start_index = 0
        while True:
            try:
                match_index = node_list[start_index:].index(search_item)
                indices_list.append(start_index + match_index)
                start_index = indices_list[-1] + 1
            except ValueError:
                # print(f'Matching indices: {indices_list}')
                return indices_list

    def grow_tree(self):
        """Extending the edit tree by a single step along the guide RNA."""

        if len(self.edit_nodes_all) == 1:
            self.root.gen_deledit_child()
            self.root.gen_nondel_children()
            self.edit_nodes_all += self.root.children
            current_index_nodes = [node for node in self.root.children]
        else:
            current_index_nodes = [node for node in self.edit_nodes_next]

        current_index_sequences = [node.sequence.seq for node in current_index_nodes]

        # if the edit node has the action 'delete' or is the root node then its children will have the same mRNA
        # index. Therefore, these will need to be generated and included in the same round of mfe determination
        # and probability calculation
        del_nodes = [node for node in current_index_nodes if node.action == 'd']
        del_sequences = [node.sequence.seq for node in del_nodes]

        # NEWER METHOD
        same_idx_c_nodes = []
        while del_nodes:
            del_nodes_unchecked = del_nodes
            del_sequences_unchecked = del_sequences

            del_nodes_checked = []
            while del_nodes_unchecked:
                test_node = del_nodes_unchecked.pop()
                test_sequence = del_sequences_unchecked.pop()

                matching_indices = self.get_match_indices(node_list=del_sequences_unchecked, search_item=test_sequence)
                nodes_to_merge = []
                for idx in reversed(matching_indices):
                    if del_nodes_unchecked[idx].action + test_node.action in ('dd', 'ip', 'pi', 'ii', 'pp'):
                        del del_sequences_unchecked[idx]
                        confirmed_match = del_nodes_unchecked.pop(idx)
                        nodes_to_merge.append(confirmed_match)
                        confirmed_match.to_merge = True
                        confirmed_match.merge_siblings = nodes_to_merge

                if nodes_to_merge:
                    test_node.to_merge = True
                    test_node.merge_siblings = nodes_to_merge
                    nodes_to_merge.append(test_node)
                    del_nodes_checked += nodes_to_merge
                else:
                    del_nodes_checked.append(test_node)

            same_idx_d_nodes = []
            for node in del_nodes_checked:
                if node.node_type is NodeType.ACTIVE:
                    node.gen_deledit_child()
                    node.gen_nondel_children()
                    same_idx_d_nodes += node.child_deledit
                    same_idx_c_nodes += node.children
                    if node.to_merge:
                        for merge_sibling in node.merge_siblings:
                            merge_sibling.node_type = NodeType.MERGED
                            merge_sibling.children = node.children
                            merge_sibling.child_deledit = node.child_deledit
                            merge_sibling.children_nondel = node.children_nondel

            del_nodes = same_idx_d_nodes
            del_sequences = [node.sequence.seq for node in del_nodes]

        c_nodes_unchecked = [node for node in same_idx_c_nodes if node.node_type is not NodeType.MERGED]
        c_sequences_unchecked = [node.sequence.seq for node in c_nodes_unchecked]

        current_index_nodes += same_idx_c_nodes
        current_index_sequences += [node.sequence.seq for node in same_idx_c_nodes]

        while c_nodes_unchecked:
            test_node = c_nodes_unchecked.pop()
            test_sequence = c_sequences_unchecked.pop()

            matching_indices = self.get_match_indices(node_list=current_index_sequences, search_item=test_sequence)
            nodes_to_merge = []
            previously_identified = []
            already_merged = False
            for idx in reversed(matching_indices):
                if current_index_nodes[idx].action + test_node.action in ('dd', 'ip', 'pi', 'ii', 'pp'):
                    confirmed_match = current_index_nodes[idx]
                    if confirmed_match.to_merge:
                        previously_identified.append(confirmed_match)
                        if confirmed_match.node_type is NodeType.MERGED:
                            already_merged = True
                    else:
                        nodes_to_merge.append(confirmed_match)
                        confirmed_match.to_merge = True

            nodes_to_merge += previously_identified

            if len(nodes_to_merge) > 1:
                for node in nodes_to_merge:
                    node.merge_siblings = nodes_to_merge
                    node.node_type = NodeType.MERGED

        if not bulk_cofold:
            self.calc_probabilities(current_index_nodes)

        next_idx_c_nodes = []
        for ci_node in current_index_nodes:
            if ci_node.node_type is NodeType.MERGED:
                if ci_node.gIndex + 1 == ci_node.guide.length:
                    ci_node.node_type = NodeType.COMPLETE
                else:
                    if not ci_node.children:
                        ci_node.gen_deledit_child()
                        ci_node.gen_nondel_children()
                        next_idx_c_nodes += ci_node.children
                    for merge_sibling in ci_node.merge_siblings:
                        if not merge_sibling.children:
                            merge_sibling.children = ci_node.children
                            merge_sibling.child_deledit = ci_node.child_deledit
                            merge_sibling.children_nondel = ci_node.children_nondel
            elif ci_node.node_type is NodeType.ACTIVE:
                if not ci_node.children:
                    ci_node.gen_deledit_child()
                    ci_node.gen_nondel_children()
                    next_idx_c_nodes += ci_node.children

        self.edit_nodes_all += same_idx_c_nodes
        self.edit_nodes_all += next_idx_c_nodes

        self.edit_nodes_next = next_idx_c_nodes



        # # NEW METHOD
        #
        # # container for all child nodes produced below
        # same_index_child_nodes = []
        #
        # # produces all the child nodes that will hold the same mIndex as the 'current_index_nodes'
        # while del_nodes:
        #     same_index_del_nodes = []
        #
        #     for edit_node in del_nodes:
        #         if (edit_node.node_type is NodeType.ACTIVE) or (edit_node.node_type is NodeType.ROOT):
        #             # create the child with 'delete' action if appropriate (based on the mRNA base)
        #             edit_node.gen_deledit_child()
        #             # append any newly created 'delete' child node - to later be iterated over for the same reason
        #             same_index_del_nodes += edit_node.child_deledit
        #             # produces remaining child nodes - 'pass' and 'insert' where appropriate
        #             edit_node.gen_nondel_children()
        #             # appends all children to a temporary list of same index child nodes. These will be added to the
        #             # main list of 'current_index_nodes' if they don't require merging with a previous node
        #             same_index_child_nodes += edit_node.children
        #
        #         # print('sleeping...1')
        #         # time.sleep(10)
        #
        #     # sets 'del_nodes' to the list of newly generated 'delete' nodes to repeat the process until there are
        #     # no same-mindex edit nodes that need creating
        #     del_nodes = same_index_del_nodes
        #
        # # container for any leaf nodes generated. These will not be merged but are tracked and later added to the
        # # 'edit_nodes_all' list for tracking all nodes produced
        # # Could also be added to the probability calculation, though currently excluded from this
        # leaf_child_nodes = []
        #
        # # check if the newly generated child nodes have the same sequence as any existing node - either in the
        # # 'current_index_nodes' list or in the list of newly generated children
        #
        # same_index_child_sequences = [node.sequence.seq for node in same_index_child_nodes]
        # first_round_checked_nodes = []
        #
        # # first check the newly generated child nodes against each other for any potential candidates for merging
        # while same_index_child_nodes:
        #     current_node = same_index_child_nodes.pop()
        #     current_seq = same_index_child_sequences.pop()
        #     if current_node.node_type is NodeType.LEAF:
        #         leaf_child_nodes.append(current_node)
        #         continue
        #     try:
        #         match_index = same_index_child_sequences.index(current_seq)
        #         match_node_action = same_index_child_nodes[match_index].action
        #
        #         # print('Found matching sequence in ChildNodes.')
        #         if (match_node_action + current_node.action[0])[-2:] in ('dd', 'ip', 'pi', 'ii', 'pp'):
        #             matching_node = same_index_child_nodes.pop(match_index)
        #             matching_sequence = same_index_child_sequences.pop(match_index)
        #
        #             merge_node = self.merge_nodes(matching_node, current_node)
        #
        #             same_index_child_nodes.append(merge_node)
        #             same_index_child_sequences.append(matching_sequence)
        #             # print('Match confirmed. Nodes merged.')
        #             # first_round_checked_nodes.append(current_node) # TEMPORARY UNTIL MERGE FUNC ACTIVATED
        #         else:
        #             # print('Actions not compatible. No merge performed.')
        #             first_round_checked_nodes.append(current_node)
        #     except ValueError:
        #         # print('No matching sequence in ChildNodes')
        #         first_round_checked_nodes.append(current_node)
        #
        # # check all child nodes - post any initial merges between themselves - against the original set of nodes at
        # # this mRNA index ('current_index_nodes'). Perform any further merges between these and append the new nodes
        # # to the list of fully checked, unique nodes
        # first_round_checked_sequences = [node.sequence.seq for node in first_round_checked_nodes]
        # fully_checked_nodes = []
        # while first_round_checked_nodes:
        #     current_node = first_round_checked_nodes.pop()
        #     current_seq = first_round_checked_sequences.pop()
        #     if current_node.node_type is NodeType.LEAF:
        #         leaf_child_nodes.append(current_node)
        #         continue
        #     try:
        #         match_index = current_index_sequences.index(current_seq)
        #         match_node_action = current_index_nodes[match_index].action
        #
        #         # print('Found matching sequence in ParentNodes.')
        #         if (match_node_action + current_node.action[0])[-2:] in ('dd', 'ip', 'pi', 'ii', 'pp'):
        #             matching_node = current_index_nodes.pop(match_index)
        #             matching_sequence = current_index_sequences.pop(match_index)
        #
        #             merge_node = self.merge_nodes(matching_node, current_node)
        #             first_round_checked_nodes.append(merge_node)
        #             first_round_checked_sequences.append(matching_sequence)
        #             # DELETE MATCHING NODE FROM ALL LISTS IT APPEARS IN
        #                 # shouldn't be necessary as nodes are merged into matching node
        #             # print('Match confirmed. Nodes merged.')
        #             # fully_checked_nodes.append(current_node) # TEMPORARY UNTIL MERGE FUNC ACTIVATED
        #         else:
        #             # print('Actions not compatible. No merge performed.')
        #             fully_checked_nodes.append(current_node)
        #     except ValueError:
        #         # print('No matching sequence in ParentNodes')
        #         fully_checked_nodes.append(current_node)
        #
        # # adds the newly created child nodes to the list of all nodes for this edit_tree
        # self.edit_nodes_all += leaf_child_nodes
        # self.edit_nodes_all += fully_checked_nodes
        #
        # # add the fully checked, child nodes to the original 'current_index_nodes' list prior to calculating mfe and
        # # probability for the current mRNA index - if 'bulk_cofold' setting is False
        # current_index_nodes += fully_checked_nodes
        #
        # if not bulk_cofold:
        #     self.calc_probabilities(current_index_nodes)
        #     # self.calc_probabilities(leaf_child_nodes)
        #
        # # generate children for all current index nodes that are ACTIVE and haven't yet had children instantiated
        # next_index_nodes = []
        # for ci_node in current_index_nodes:
        #     if ci_node.node_type is NodeType.ACTIVE:
        #         if not ci_node.children:
        #             ci_node.gen_deledit_child()
        #             ci_node.gen_nondel_children()
        #             next_index_nodes += ci_node.children
        #
        # # these 'next_index_nodes' are appended to 'edit_nodes_all' and the ACTIVE nodes among them are set to the
        # # 'edit_nodes_next' list, for use in the next round of 'edit_tree' growth
        # self.edit_nodes_all += next_index_nodes
        # self.edit_nodes_next = next_index_nodes



        # # OLD METHOD
        #
        # next_nodes = current_index_nodes[:]
        #
        # while next_nodes:
        #     new_deledit_nodes = []
        #     new_deledit_sequences = []
        #
        #     for edit_node in next_nodes:
        #         if (edit_node.node_type is NodeType.ACTIVE) or (edit_node.node_type is NodeType.ROOT):
        #             edit_node.gen_deledit_child()
        #             if edit_node.child_deledit:
        #                 try:
        #                     orig_index = current_index_sequences.index(edit_node.child_deledit[0].sequence.seq)
        #                     orig_node = current_index_nodes[orig_index]
        #                     print('Found matching sequence')
        #                     if orig_node.action == edit_node.action:
        #                         self.merge_nodes(orig_node, edit_node)
        #                     else:
        #                         print('Not compatible for merging')
        #                 except ValueError:
        #                     print('No sequences to merge')
        #                 new_deledit_nodes.append(edit_node.child_deledit[0])
        #                 new_deledit_sequences.append(edit_node.child_deledit[0].sequence.seq)
        #
        #     self.edit_nodes_all += new_deledit_nodes
        #     next_nodes = new_deledit_nodes
        #     current_index_nodes += new_deledit_nodes
        #     current_index_sequences += new_deledit_sequences
        #
        # self.edit_nodes_next = current_index_nodes
        # if not bulk_cofold:
        #     self.calc_probabilities(self.edit_nodes_next)
        #
        # next_index_nodes = []
        # for edit_node in current_index_nodes:
        #     if (edit_node.node_type is NodeType.ACTIVE) or (edit_node.node_type is NodeType.ROOT):
        #         edit_node.gen_nondel_children()
        #         next_index_nodes += edit_node.children_nondel
        #
        # self.edit_nodes_all += next_index_nodes
        # self.edit_nodes_next = next_index_nodes

        # BELOW IS CONSISTENT FOR ALL METHODS

        # the edit levels count is increased by 1, corresponding to the number of 'edit_tree' growth rounds completed
        self.edit_levels += 1

        # if there are no newly produced nodes to be appended to 'edit_nodes_next' then the 'edit_tree' has reached
        # its terminal state and no further growth is possible. It is therefore set to complete
        if not self.edit_nodes_next:
            self.is_complete = True
            # print(f'Edit tree completed with {self.edit_levels} levels, {len(self.edit_nodes_all)} nodes.')

    def select_progressed_nodes(self) -> list:
        """Takes all the COMPLETE nodes from the edit tree and selects the sequences to progress to the subsequent
        round of guide RNA docking.
        Selection is based on the probability product of the node, its calculated MFE, and the number of
        sequences to progress, as defined in the run settings.

        Full order of ranking criteria is as follows:
            - probability_product (descending)
            - mfe (ascending)
            - mismatches (ascending)
            - gIndex (descending)
            - edit_level (ascending)"""

        # select the edit_nodes which are classified as 'COMPLETE'. These edit_nodes have completed the editing
        # process, reaching the end of the guide RNA template without becoming a 'LEAF' node
        complete_nodes = [edit_node for edit_node in self.edit_nodes_all if edit_node.node_type is NodeType.COMPLETE]

        # if no such 'COMPLETE' nodes exist, an empty list is returned
        if not complete_nodes:
            return []

        # the number of sequences to progress to the subsequent round of guide RNA docking is no bigger than the
        # upper limit, set out as 'sequences_to_progress' in the 'run_settings' module
        number_to_progress = min(sequences_to_progress, len(complete_nodes))

        # create a list of dicts, where each dict contains the full list of variables for each complete_node
        comp_vars = [vars(node) for node in complete_nodes]
        # the list of dicts is used to make a DataFrame, enabling simple, multi-parameter sorting
        comp_df = pd.DataFrame(comp_vars)
        comp_df = comp_df[comp_df['mfe'] < min_mfe_to_progress]

        ranked = comp_df.sort_values(by=['probability_product', 'mfe', 'mismatches', 'gIndex', 'edit_level'],
                                     ascending=[False, True, True, False, True])
        # original index is replaced by the ranked index. The ranked index is used to select the top n edit_nodes.
        # The original index of these selected edit nodes can then be used to retrieve the edit_nodes from the
        # 'complete_nodes' list, since its index matches the original DataFrame index
        candidate_indices = ranked.reset_index(drop=False).iloc[:number_to_progress]['index']

        # chosen edit_nodes are selected and placed in a new 'progressed_nodes' list, which is then returned
        progressed_nodes = []
        for index in candidate_indices:
            selected_node = complete_nodes[index]
            # chosen 'edit_node' is updated to show that it was progressed to the next round of guide_node creation
            selected_node.progressed = True
            progressed_nodes.append(selected_node)

        return progressed_nodes
