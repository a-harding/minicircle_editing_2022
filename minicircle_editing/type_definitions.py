"""Definition of the various type classes used for specifying settings of a program run and some node/sequence
attributes."""


from enum import Enum, auto


class SequenceType(Enum):
    """Defines a sequence as a guide or messenger RNA."""
    MESSENGER = auto()
    GUIDE = auto()


class DockingMode(Enum):
    """How the gRNA selection may be weighted in favour of those docking at specific sites."""
    INITIATION = auto()  # Initial gRNA selection weighted in favour of gRNAs at 3' end of mRNA
    CURRENT_SITE = auto()  # Subsequent gRNA selection favours docking near previous edit site
    INITIATION_AND_CURRENT = auto()  # Both of the above
    NO_WEIGHTING = auto()


class gRNAExclusion(Enum):
    """Determines if gRNAs that have previously edited the transcript can be selected for further editing."""
    ALL = auto()  # All previously used gRNAs are excluded from further binding
    ONE = auto()  # The most recently used gRNA is excluded
    NONE = auto()  # No gRNAs are excluded. All are considered without preference


class EditMode(Enum):
    """Sets the type of editing that will be used for the model."""
    STEP = auto()  # Edit nodes are created in the manner of a binary decision tree
                        # Co-folding is performed and probabilities calculated after generating each pool of child nodes
                        # Progression threshold set by probability and total editing errors
    BULK = auto()  # Edit nodes are generated in similar binary manner but no co-folding is performed until
                        # all possible edit nodes have been produced for the entire gRNA - sent as a batch
                        # Threshold on total editing errors to prevent needless node generation


class CofoldMode(Enum):
    """Setting how much of a guide sequence to include when sending to RNAcofold. mRNA sequence length matches
    that of the guide RNA."""
    WHOLE_GUIDE = auto() # Includes the entire guide RNA sequence
    TO_INDEX = auto() # The 5' end of the guide, up to the current guide base being considered in editing
    TO_INDEX_PLUS = auto() # As above, plus an additional number of bases the size of the editing window
    WINDOW_CENTRED = auto() # A region the size of the editing window, centred on the current gRNA editing base


class NodeType(Enum):
    """The different classifications of an edit node."""
    ROOT = auto() # The first node in a new edit tree
    ACTIVE = auto() # A non-terminal node
    LEAF = auto() # A terminal node. It has exceeded a limiting parameter
    COMPLETE = auto() # A node which has reached the end of the guide RNA and its editing is complete
    MERGED = auto() # Node has been used to form a merged node with one or more other edit nodes. No longer considered
    # for further operations in the edit tree, but not necessarily a LEAF


class OutputNodes(Enum):
    """Selection of nodes to include in the output CSV files."""
    ALL = auto() # All edit nodes of all guide nodes.
    COMPLETE = auto() # Edit nodes that completed editing.
    PROGRESSED = auto() # Only nodes that progressed to the next round of editing.
