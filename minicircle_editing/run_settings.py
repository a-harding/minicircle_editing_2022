"""User-specified settings for an individual run. Primary means of controlling run-time and output of the model."""


from type_definitions import DockingMode, gRNAExclusion, EditMode, CofoldMode


# edit_mode = EditMode.STEP  # which mode of editing will be used for node generation
no_of_grnas_first = 1  # how many gRNAs to investigate with unedited sequence
min_no_grnas_subsequent = 2  # how many gRNAs to investigate with partially edited sequences
max_no_grnas_subsequent = 5
guides_to_cofold = 50 # number of guides to send to RNAcofold during gRNA selection step
cofold_mode = CofoldMode.TO_INDEX_PLUS # how much of the gRNA to send to RNAcofold
guide_end_allowance = 3 # how many bases to disregard when calculating mismatches near the terminal end of a guide
sequences_to_progress = 1  # how many edited sequences to take forward to next gRNA editing
mismatch_threshold_anchor = 2  # no. of mismatches allowed when identifying gRNA anchors
mismatch_threshold_editing = 2  # no. of mismatches allowed during editing before node becomes leaf
editing_window = 7  # no. of mRNA bases beyond current bp to include in RNAcofold
max_anchor = 15  # maximum permissible anchor length when performing docking
min_anchor = 8  # minimum permissible anchor length when docking
compare_all_nodes = True  # compare all nodes of the same level when calculating probability?
probability_threshold = 0.01  # if above bool is False, set to 0.1. If True, set lower - 0.01, 0.001
docking_mode = DockingMode.CURRENT_SITE  # gRNAs preferentially dock to 3' end, near prev. site, or no pref.
previous_gRNA_exclusion = gRNAExclusion.ALL # how gRNAs previously used to edit the sequence should be considered in docking
bulk_cofold = True # If True, edit trees are built with mismatches as the only threshold, with cofolding performed
                    # on all complete, non-leaf nodes once the tree is fully built
maximum_edit_graph_nodes = 1000
short_sequence_editing = True
proportion_to_dock = 0.5
minimum_mfe = -7
minimum_adjusted_mfe = -7
min_mfe_to_progress = -30