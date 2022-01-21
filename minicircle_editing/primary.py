""""""


import time
from sequence_import import import_seqeunces
from log_gen import gen_setup_file
from docking import select_guides
from guide_trees import GuideTree
from notifications import email_notification

def main():
    """Main container for executing the entire model."""

    messenger, guides_dict, gene, paths, edited_seq = import_seqeunces() # imports all initial sequences
    start = time.perf_counter()

    log_path = gen_setup_file(gene, paths) # save a log file with chosen settings for the current run

    # perform initial docking and select starting sequences
    initial_duplexes = select_guides(messenger=messenger, guides_dict=guides_dict, initial=True)
    # each duplex is list with order [guide_name, mRNA_dock, mRNA_sequence, guide_index]
    time.sleep(2)

    for duplex in initial_duplexes:
        GuideTree(initial_duplex=duplex, guides_dict=guides_dict, gene=gene, log_path=log_path, edited_seq=edited_seq)



    # email_notification('Hi')