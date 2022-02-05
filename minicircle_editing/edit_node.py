"""Defines the class of edit node and associated functions. Each node is an editing decision at a single base."""


from sequence_import import Sequence
from docking import gen_cofold_string
from type_definitions import NodeType
from run_settings import mismatch_threshold_editing, probability_threshold, guide_end_allowance
import RNA
from copy import copy


class EditNode:
    """Class representing a single editing decision during the editing process."""

    mismatch_set = {'gg', 'cc', 'aa', 'uu', 'ga', 'ag', 'ac', 'ca', 'cu', 'uc'}

    def __init__(self, edit_tree, action: str, parent=None):
        self.edit_tree = edit_tree
        self.guide = self.edit_tree.guide_sequence
        self.parent = parent
        self.parents = [self.parent]
        self.merge_parents = []
        self.merge_siblings = []
        self.children = []
        self.child_deledit = []
        self.children_nondel = []
        self.action = action
        self.progressed = False
        self.to_merge = False
        self.probability_calculated = False

        self.edit_level: int
        self.action_log: str
        self.init_sequence: Sequence
        self.sequence: Sequence
        self.mIndex: int
        self.gIndex: int
        self.node_type: NodeType
        self.probability: float
        self.probability_product: float
        self.combined_parent_prob_prod: float
        self.cofold_string: str
        self.mfe: float
        self.alignment: str

        # self.pairs: str

    def __lt__(self, other):
        return self.edit_level < other.edit_level

    def get_next_index_step(self):
        """Checks if the next base to be considered for editing is the same as the current index or the current
        index + 1, based on the current node's editing action."""

        # first check if the next mIndex is advanced from current mIndex. If the current node's action was to 'pass'
        # or 'insert' at the current base, then the next base considered will be current mIndex + 1
        if self.action in {'p', 'i'}:
            next_index_step = 1
        # if the current node's action was to 'delete' or perform no action (because it is a root node), then the
        # next base considered will be at the current mIndex
        else:
            next_index_step = 0

        return next_index_step

    def gen_nondel_children(self) -> None:
        """Instantiates remaining child nodes if they have yet to be produced. Stores all in self.children."""

        # if the node has already generated its child nodes, no further children are created
        if self.children:
            pass
        else:
            # first check if the next mIndex is advanced from current mIndex.
            next_index_step = self.get_next_index_step()

            # all non-leaf nodes produce a 'pass' node, in which no editing of the mRNA sequence will occur and
            # the only action will be to advance the editing machinery by a single base
            pass_node = EditNodeChild(edit_tree=self.edit_tree, action='p', parent=self)
            self.children_nondel.append(pass_node)

            # if the most recent action was to delete a 'U', no insert will be performed
            if self.action != 'd':
                # a further check (indirectly checking that no sibling has been created with the 'delete' action) to
                # confirm that the current base isn't a 'U'. Otherwise no 'insert' child will be created
                if self.sequence.seq[self.mIndex + next_index_step] != 'u':
                    insert_node = EditNodeChild(edit_tree=self.edit_tree, action='i', parent=self)
                    self.children_nondel.append(insert_node)

        # once all children are created in their sub-categories, they are combined into the self.children list
        self.children = self.children_nondel + self.child_deledit

        return None

    def gen_deledit_child(self) -> None:
        """Instantiates the child node with a 'delete' action if appropriate."""

        # first check if the next mIndex is advanced from current mIndex.
        next_index_step = self.get_next_index_step()

        # if the next base considered is a 'u' then one edit option is to delete that 'u'. In which case a child node
        # with the 'delete' action is created
        if self.action != 'i':
            if self.sequence.seq[self.mIndex + next_index_step] == 'u':
                new_node = EditNodeChild(edit_tree=self.edit_tree, action='d', parent=self)
                self.child_deledit.append(new_node)
            else:
                pass

        return None

    def count_all_mismatches(self) -> int:
        """Counts the number of mismatches in the entire sequence paired between the guide RNA and mRNA.

        A slower method that iterates over the entire paired sequence. Only used to confirm the correct
        functionality of the cumulative use of count_mismatch method. Not used in the final program."""

        pair_string = ''
        mismatches = 0
        for i in range(self.gIndex + 1):
            m_base = f'{self.sequence.seq[i + self.edit_tree.init_mIndex - self.edit_tree.init_gIndex]}'
            g_base = f'{self.guide.seq[i]}'
            current_pair = m_base + g_base
            if current_pair in EditNode.mismatch_set:
                mismatches += 1
            pair_string = pair_string + f'({current_pair})'

        # self.pairs = pair_string

        return mismatches

    def set_cofold_string(self) -> str:
        """Calls the gen_cofold_string function from the docking module and returns the correctly formatted string
        for sending to RNAcofold."""

        cofold_string = gen_cofold_string(messenger_sequence=self.sequence.seq,
                                          guide_sequence=self.guide.seq,
                                          docking_idx=self.edit_tree.init_mIndex - self.edit_tree.init_gIndex,
                                          gIndex=self.gIndex)

        return cofold_string

    def cofold_sequence(self) -> float:
        """Calculates MFE by sending the correctly formatted cofold string to RNAcofold."""

        self.alignment, self.mfe = RNA.cofold(self.cofold_string)

        return self.mfe

    def set_node_type(self) -> NodeType:
        """Sets the initial state of the node type based on mismatches and index positions."""

        if not self.parent:
            return NodeType.ROOT
        elif self.mismatches > mismatch_threshold_editing:
            if self.gIndex + 1 + guide_end_allowance >= self.guide.length:
                return NodeType.COMPLETE
            else:
                return NodeType.LEAF
        elif (self.gIndex + 1 == self.guide.length) or (self.mIndex + 1 == self.sequence.length):
            return NodeType.COMPLETE
        else:
            return NodeType.ACTIVE

    def reset_node_type(self):
        """Checks the thresholds to determine if the node is terminal or active."""

        if self.node_type is NodeType.ACTIVE:
            if self.probability_product < probability_threshold:
                self.node_type = NodeType.LEAF


class EditNodeRoot(EditNode):
    """Subclass of EditNode. Called when producing the initial node for a new edit tree. Sets the starting
    probabilities, sequences, and indices."""

    def __init__(self, edit_tree, action: str, parent=None, init_seq=None):
        super().__init__(edit_tree=edit_tree, action=action)

        self.id = f'{self.edit_tree.id}_{self.action}'
        self.parent_id = None
        self.edit_level = 0
        self.action_log = self.action

        self.mIndex = copy(self.edit_tree.init_mIndex)
        self.gIndex = copy(self.edit_tree.init_gIndex)

        self.init_sequence = init_seq
        self.sequence = self.init_sequence

        self.prev_gIndex_mismatch_total = 0
        self.mismatches = 0

        self.probability = 1.0
        self.probability_product = 1.0

        self.node_type = self.set_node_type()
        self.cofold_string = self.set_cofold_string()
        # self.mismatches_full = self.count_all_mismatches()


class EditNodeChild(EditNode):
    """Used for all edit nodes in an edit tree that aren't the first node. Carries all the editing logic and contain
    the progressively edited mRNA sequences."""

    def __init__(self, edit_tree, action: str, parent: EditNode):
        super().__init__(edit_tree=edit_tree, action=action, parent=parent)

        self.id = f'{self.parent.id}{self.action}'
        self.parent_id = self.parent.id
        self.edit_level = self.parent.edit_level + 1
        self.action_log = self.parent.action_log + self.action

        if self.parent.action in {'p', 'i'}:
            self.mIndex = self.parent.mIndex + 1
            self.gIndex = self.parent.gIndex + 1
        else:
            self.mIndex = copy(self.parent.mIndex)
            self.gIndex = copy(self.parent.gIndex)
        # print(f'Action: {action}; mIndex - Parent: {self.parent.mIndex}, Self: {self.mIndex}')

        self.init_sequence = self.parent.sequence
        self.sequence = self.edit(init_sequence=self.init_sequence)

        if self.parent.action == 'd':
            self.prev_gIndex_mismatch_total = self.parent.prev_gIndex_mismatch_total
        else:
            self.prev_gIndex_mismatch_total = self.parent.mismatches
        self.mismatches = self.prev_gIndex_mismatch_total + self.check_mismatch()

        self.node_type = self.set_node_type()
        self.cofold_string = self.set_cofold_string()
        # self.mismatches_full = self.count_all_mismatches()

    def edit(self, init_sequence: Sequence) -> Sequence:
        """Performs the editing decision on the inherited RNA sequence."""

        prev_seq = init_sequence.seq3to5
        if self.action == 'i':
            new_seq = f'{prev_seq[:self.mIndex]}U{prev_seq[self.mIndex:]}'
        elif self.action == 'd':
            new_seq = f'{prev_seq[:self.mIndex]}{prev_seq[self.mIndex + 1:]}'
        else:
            new_seq = prev_seq

        return Sequence(name=self.id, sequence=new_seq, seq_is_5to3=False)

    def check_mismatch(self) -> int:
        """Checks if the current base pair is a mismatch."""

        current_base_m = self.sequence.seq3to5[self.mIndex]
        current_base_g = self.edit_tree.guide_sequence.seq[self.gIndex]

        if f'{current_base_m}{current_base_g}' in EditNode.mismatch_set:
            return 1
        else:
            return 0
