"""This module is responsible for selection of the fasta files to import sequences from and for the
act of importing all sequences (both mRNA and gRNA) to be used in the model.

It also defines a simple Sequence class for containing the important features of these sequences."""


from pathlib import Path
from os import listdir
from type_definitions import SequenceType


class Sequence:
    """Simple sequence class for containing information about each mRNA or gRNA."""

    def __init__(self, name: str, sequence: str, seq_is_5to3=True, seq_type=SequenceType.MESSENGER):
        self.seq_type = seq_type
        self.name = name
        if seq_is_5to3:
            self.seq5to3 = sequence.lower()
        else:
            self.seq5to3 = sequence[::-1].lower()
        self.seq3to5 = self.seq5to3[::-1]
        self.length = len(sequence)
        self.seq: str # Standard sequence for comparisons. Set to 3' to 5' if mRNA, or 5' to 3' for gRNA
        if self.seq_type is SequenceType.MESSENGER:
            self.seq = self.seq3to5
        else:
            self.seq = self.seq5to3

    # def __eq__(self, other):
    #     return self.seq5to3 == other.seq5to3


def import_mrna(filepath: Path) -> Sequence:
    """Imports the mRNA name and sequence from a fasta file. This is then assigned to a Sequence object.
    Returns a single mRNA Sequence object."""

    with open(filepath) as file:
        contents = file.read()
        whole_file = contents.strip().lower().split('\n')
        mrna_name = whole_file[0]
        mrna_sequence = whole_file[1]

        # Place reference mRNA into a Sequence object
        messenger_sequence = Sequence(mrna_name, mrna_sequence, seq_is_5to3=True, seq_type=SequenceType.MESSENGER)

    return messenger_sequence


def import_grnas(filepath: Path) -> dict:
    """Imports all gRNAs from a fasta file, splitting them out into separate Sequence objects.
    Return a list of gRNA Sequence objects."""

    with open(filepath) as file:
        contents = file.read()
        all_grnas = contents.strip().lower().split('\n')

    # Split out the name data and sequence data of gRNAs into separate lists
    grna_names = all_grnas[0::2]
    grna_sequences = all_grnas[1::2]

    # Make a list of all gRNA Sequence objects
    guides = [Sequence(name, seq, seq_type=SequenceType.GUIDE) for name, seq in zip(grna_names, grna_sequences)]
    guides_dict = dict(zip(grna_names, guides))

    return guides_dict


def initialise_sequences(filepath_mRNA: Path, filepath_gRNA: Path, filepath_edited_mRNA: Path) -> tuple:
    """Calls to both mRNA and gRNA importers to establish all working sequences."""

    unedited_sequence = import_mrna(filepath_mRNA)
    guides_dict = import_grnas(filepath_gRNA)
    edited_sequence = import_mrna(filepath_edited_mRNA)

    return unedited_sequence, guides_dict, edited_sequence


def user_selection(file_names: list) -> str:
    """User selects option (gene and associated gRNA file) from presented lists."""

    selection = len(file_names)
    for i in range(len(file_names)):
        print(f'{i}. {file_names[i]}')
    while True:
        try:
            selection = int(input())
        except ValueError:
            print('Please enter an integer.')
        if (selection < len(file_names)) and (selection >= 0):
            break
        else:
            print('Selection not in range. Choose again:')
    print(f'Selected: {file_names[selection]}')

    return file_names[selection]


def get_files(mrna_folder: Path, grna_folder: Path, genes: list, edited_mrna_folder: Path) -> tuple:
    """User selection of the gene to investigate. Returns the appropriate mRNA file name."""

    mrna_file_names = sorted(listdir(mrna_folder))
    grna_file_names = sorted(listdir(grna_folder))
    edited_mrna_file_names = sorted(listdir(edited_mrna_folder))

    print('Select the unedited mRNA sequence file you wish to begin editing:')

    gene_file_name = user_selection(mrna_file_names)
    output_filepath_mrna = mrna_folder / gene_file_name

    selected_gene = ''
    for gene in genes:
        if gene.lower() in gene_file_name.lower():
            selected_gene = gene

    potential_grna_files = []
    for file in grna_file_names:
        if selected_gene.lower() in file.lower():
            potential_grna_files.append(file)

    if len(potential_grna_files) > 1:
        output_filepath_grna = grna_folder / user_selection(potential_grna_files)
    else:
        output_filepath_grna = grna_folder / potential_grna_files[0]
        print(f'Selected: {potential_grna_files[0]}')

    potential_edited_mrna_files = []
    for file in edited_mrna_file_names:
        if selected_gene.lower() in file.lower():
            potential_edited_mrna_files.append(file)

    if len(potential_edited_mrna_files) > 1:
        output_filepath_edited_mrna = edited_mrna_folder / user_selection(potential_edited_mrna_files)
    else:
        output_filepath_edited_mrna = edited_mrna_folder / potential_edited_mrna_files[0]
        print(f'Selected: {potential_edited_mrna_files[0]}')

    return output_filepath_mrna, output_filepath_grna, selected_gene, output_filepath_edited_mrna


def define_import_filepaths() -> tuple:
    """User-assisted selection of the filepaths for sequence import.
    Filepaths are interpreted from gene name selection."""

    genes = ['A6', 'COX3', 'CR3', 'CR4', 'CYB', 'MURF2', 'ND3', 'ND7', 'ND8', 'ND9', 'RPS12']

    folderpath_mrna = Path('Sequences/mRNAs_unedited/')
    folderpath_grna = Path('Sequences/gRNAs/')
    folderpath_edited_mrna = Path('Sequences/mRNAs_edited')

    filepath_mrna, filepath_grna, selected_gene, filepath_edited_mrna = get_files(
        mrna_folder=folderpath_mrna, grna_folder=folderpath_grna, genes=genes, edited_mrna_folder=folderpath_edited_mrna
    )

    return filepath_mrna, filepath_grna, selected_gene, filepath_edited_mrna


def import_seqeunces() -> tuple:
    """Wrapper for calling functions selecting sequence files, then importing sequences from the
    selected locations. Returns the reference mRNA and guide RNA Sequence objects"""

    filepath_mrna, filepath_guides, selected_gene, filepath_edited_mrna = define_import_filepaths()
    messenger_sequence, guides_dict, edited_sequence = initialise_sequences(filepath_mRNA=filepath_mrna,
                                                                            filepath_gRNA=filepath_guides,
                                                                            filepath_edited_mRNA=filepath_edited_mrna)

    return messenger_sequence, guides_dict, selected_gene, [filepath_mrna, filepath_guides], edited_sequence
