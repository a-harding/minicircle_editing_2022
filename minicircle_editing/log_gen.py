"""Responsible for generating, saving, and updating run log files."""


from pathlib import Path
from datetime import date
import time
from run_settings import *


def gen_setup_file(gene: str, paths: list[Path]) -> Path:
    """Produces a .txt file recording the initialisation features of the run."""

    start_date = date.today()
    start_time = time.strftime("%H-%M-%S", time.localtime())

    dest_dir = Path.home() / Path(f'Documents/HPGH/SchnauferLab/Data/{gene}/{start_date}_{start_time}')
    dest_dir.mkdir(parents=True, exist_ok=True)

    log_filename = f'initiator_settings.txt'
    log_path = dest_dir / Path(f'{log_filename}')

    path_m, path_g = paths

    with open(log_path, mode='w', encoding='utf-8') as f:
        f.write(f'*****\nGENE AND FILE SELECTION\n*****\n'
                f'Gene: {gene}\n'
                f'Path to unedited mRNA sequence: {path_m}\n'
                f'Path to gRNA pool: {path_g}\n'
                f'\n*****\ngRNA SELECTION\n*****\n'
                f'Number of gRNAs investigated (first round, subsequent rounds): {no_of_grnas_first}, '
                f'{min_grnas_subsequent}\n'
                f'Normalisation of docking by proximity to editing site: {docking_mode}\n'
                f'How many previously used gRNAs are excluded from subsequent gRNA selection: {previous_gRNA_exclusion}\n'
                f'Maximum number of sequences progressed from a single gRNA\'s editing tree: {sequences_to_progress}\n'
                f'\n*****\nANCHOR DETERMINATION\n*****\n'
                f'Maximum length of anchor: {max_anchor}\n'
                f'Mismatches allowed in anchor region: {mismatch_threshold_anchor}\n'
                f'\n*****\nEDITING\n*****\n'
                f'Mismatches allowed in editing: {mismatch_threshold_editing}\n'
                f'Probability threshold, below which a node is converted to a leaf: {probability_threshold}\n'
                f'RNAcofold mode: {cofold_mode}\n'
                f'Length of editing window: {editing_window}\n'
                f'Bulk cofold mode active: {bulk_cofold}'
                f'\n*****\nRUN LOG\n*****\n')

    return log_path


def append_log(log_path: Path, new_events: list[str]) -> None:
    """Appends a given list of events to the initial log file."""

    with open(log_path, mode='a', encoding='utf-8') as f:
        f.writelines(new_events)

    return None


def convert_timestamp(time_start: float, time_end: float) -> tuple:
    """Converts two timestamp floats generated by time.perf_counter() into hours, minutes, seconds for logging."""
    elapsed = time_end - time_start
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60

    return hours, minutes, seconds
