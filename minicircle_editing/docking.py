"""Performs the docking of guide RNAs to the mRNA sequence and selects the appropriate guide(s) to progress
to the stage of editing."""
import pprint
import time

import numpy as np
import pandas as pd
import scipy.stats as st
import RNA
from sequence_import import Sequence
from type_definitions import CofoldMode, DockingMode, gRNAExclusion
from run_settings import (
    no_of_grnas_first, no_of_grnas_subsequent, editing_window, docking_mode, max_anchor, min_anchor, guides_to_cofold,
    previous_gRNA_exclusion, cofold_mode, proportion_to_dock, minimum_mfe, guide_end_allowance,
    mismatch_threshold_anchor as mismatches_allowed
)


def split_convert(sequence: str) -> np.array:
    """Converts a sequence string to a numeric numpy array for matrix calculation of binding."""

    conversion_dict = {'G': 1, 'C': 2, 'A': 11, 'U': 12,
                       'g': 1, 'c': 2, 'a': 11, 'u': 12}
    converted = np.array([conversion_dict[char] for char in sequence])

    return np.array(converted)


def align_guide(messenger: np.array, guide: np.array) -> dict:
    """Aligns the guide to the reference mRNA sequence by projecting both sequences together into a matrix."""

    # if len(guide.shape) == 1:
    #     guide = guide.reshape(len(guide), 1)

    matrix = np.where((np.abs(messenger - guide) % 10) == 1, 0, 1)

    for i in range(1, guide.size):
        matrix[i] = np.hstack((matrix[i, i:], np.full(i, 9)))

    summed = np.cumsum(matrix, axis=0)
    guideIdx, refIdx = np.where(summed <= mismatches_allowed)

    return dict(zip(refIdx, guideIdx))


def align_all_guides(messenger_sequence: str, guides_dict: dict, excluded_guides: list) -> dict:
    """Align all guides against the given reference sequence and select the candidates for RNAcofold."""

    mes = split_convert(messenger_sequence)

    alignments_dict = {}
    for guide in guides_dict.values():
        if guide.name in excluded_guides:
            pass
        else:
            alignments_dict[guide.name] = align_guide(mes, split_convert(guide.seq[:max_anchor]).reshape(max_anchor, 1))

    alignments_df = pd.DataFrame(alignments_dict)
    nth_guide_value = np.sort(alignments_df.to_numpy(), axis=None)[::-1][guides_to_cofold - 1]

    selected_indices = {}
    for name in guides_dict.keys():
        if name in excluded_guides:
            pass
        else:
            indices = alignments_df.index[alignments_df[name] >= nth_guide_value]
            if len(indices):
                selected_indices[name] = indices

    return selected_indices


def gen_cofold_string(messenger_sequence: str, guide_sequence: str, docking_idx: int, gIndex=None) -> str:
    """Produces the string for sending to RNAcofold, with sequence trimming in line with the run settings for
    how much of the guide to include in co-folding."""

    half_window = editing_window // 2

    if gIndex:
        guide_index = min(gIndex, len(guide_sequence) - 1)
        if cofold_mode is CofoldMode.WHOLE_GUIDE:
            guide_trimmed = guide_sequence
        elif cofold_mode is CofoldMode.TO_INDEX:
            guide_trimmed = guide_sequence[:guide_index + 1]
        elif cofold_mode is CofoldMode.TO_INDEX_PLUS:
            guide_trimmed = guide_sequence[:min(guide_index + 1 + editing_window, len(guide_sequence) - 1)]
        elif cofold_mode is CofoldMode.EDITING_WINDOW:
            idx_lower = max(0, guide_index + 1 - half_window)
            idx_upper = min(guide_index + half_window, len(guide_sequence) - 1)
            guide_trimmed = guide_sequence[idx_lower:idx_upper]

            midx_lower = max(0, docking_idx + gIndex + 1 - half_window)
            midx_upper = min(docking_idx + gIndex + half_window, len(messenger_sequence) - 1)
            mrna_trimmed = messenger_sequence[midx_lower:midx_upper]

    else:
        trimmed_index = int(len(guide_sequence) * proportion_to_dock)
        guide_trimmed = guide_sequence[:trimmed_index]

    if (not gIndex) or (cofold_mode is not CofoldMode.EDITING_WINDOW):
        midx_lower = docking_idx
        midx_upper = min(len(guide_trimmed) + docking_idx, len(messenger_sequence))
        mrna_trimmed = messenger_sequence[midx_lower:midx_upper]

    mrna_trimmed = mrna_trimmed[::-1]

    # print(f'G: {guide_trimmed} - {len(guide_trimmed)}')
    # print(f'M: {mrna_trimmed} - mI {mIndex}, full_len {len(messenger_sequence)} - {len(mrna_trimmed)}')

    cofold_string = f'{mrna_trimmed}&{guide_trimmed}'
    # CHECKING IF THERE IS AN ERROR IN PRODUCING THE COFOLD STRING - BASED ON AN INDEX ERROR
    if cofold_string[0] == '&':
        print('No mRNA string!\nSleeping...')
        for i in range(60):
            print(60 - i)
            time.sleep(1)

    return cofold_string


def determine_mfe(messenger_sequence: str, guides_dict: dict, indices_dict: dict) -> pd.DataFrame:
    """Calculate MFE using RNAcofold prediction algorithm."""

    MFE_pair_list = []
    for name, indices in indices_dict.items():
        MFE_temp_list = []
        for idx in indices:
            cofold_string = gen_cofold_string(messenger_sequence, guides_dict[name].seq, idx)
            alignment, MFE = RNA.cofold(cofold_string)
            MFE_temp_list.append([name, idx, MFE])
        MFE_pair_list += MFE_temp_list

    return pd.DataFrame(MFE_pair_list, columns=['Guide_name', 'mDock', 'MFE'])


def normalise_mfe(MFE_df: pd.DataFrame, current_mIndex) -> pd.DataFrame:
    """Adjusts MFE to favour binding near to the designated mRNA index."""

    mfes = MFE_df['MFE'].to_numpy()
    idxs = MFE_df['mDock'].to_numpy()

    z_score = abs((current_mIndex - idxs) / (editing_window * 2))
    normalisation_factor = (1 - st.norm.cdf(z_score)) * 2
    adjusted_mfes = mfes * normalisation_factor

    MFE_df['Normalisation_factor'] = normalisation_factor
    MFE_df['Adjusted_MFE'] = adjusted_mfes

    return MFE_df


def sort_candidates(mfe_df: pd.DataFrame, current_mIndex=0, initial=False) -> pd.DataFrame:
    """Determines the need for generating a normalised MFE based on run settings. Then sorts the candidate guides
    by the MFE/adjusted MFE accordingly. Returns the sorted DatFrame"""

    if initial:
        if (docking_mode is DockingMode.INITIATION) or (docking_mode is DockingMode.INITIATION_AND_CURRENT):
            adjusted_mfe = normalise_mfe(mfe_df, current_mIndex)
        else:
            return mfe_df.sort_values(by=['MFE'])
    else:
        if (docking_mode is DockingMode.CURRENT_SITE) or (docking_mode is DockingMode.INITIATION_AND_CURRENT):
            adjusted_mfe = normalise_mfe(mfe_df, current_mIndex)
        else:
            return mfe_df.sort_values(by=['MFE'])

    return adjusted_mfe.sort_values(by=['Adjusted_MFE'])


def get_index(messenger_sequence: str, guide_sequence: str, mIndex: int):
    """Aligns the candidate sequences and determines the base to begin editing - the gIndex."""

    match_set = {'gc', 'cg', 'au', 'ua', 'gu', 'ug'}
    anchor_length = 0
    mismatches = 0
    consecutive_mismatches = 0

    for m, g in zip(messenger_sequence[mIndex:mIndex + len(guide_sequence)].lower(), guide_sequence.lower()):
        if m + g in match_set:
            consecutive_mismatches = 0
        else:
            mismatches += 1
            consecutive_mismatches += 1
        anchor_length += 1

        if mismatches == mismatches_allowed:
            break

    gIndex = anchor_length - consecutive_mismatches

    if gIndex + guide_end_allowance >= len(guide_sequence):
        return 0

    return gIndex


def get_excluded_guides(previous_guides: list) -> list:
    """Determines which guides to exclude based on the run settings."""

    if previous_guides:
        if previous_gRNA_exclusion is gRNAExclusion.ALL:
            return previous_guides
        elif previous_gRNA_exclusion is gRNAExclusion.ONE:
            return previous_guides[-1]
        elif previous_gRNA_exclusion is gRNAExclusion.NONE:
            return []
    else:
        return []


def select_guides(messenger: Sequence, guides_dict: dict, previous_guides=None, current_mIndex=0, initial=False):
    """Selection of the best guide RNA candidates for the given sequence."""

    if previous_guides:
        if len(guides_dict) == len(previous_guides):
            return []

    excluded_guides = get_excluded_guides(previous_guides)
    messenger_sequence = messenger.seq
    # print(messenger_sequence)
    indices = align_all_guides(messenger_sequence=messenger_sequence,
                               guides_dict=guides_dict,
                               excluded_guides=excluded_guides)
    mfe_df = determine_mfe(messenger_sequence, guides_dict, indices)
    candidates_sorted = sort_candidates(mfe_df, current_mIndex, initial=initial)

    mfe_filtered = candidates_sorted[candidates_sorted['MFE'] < minimum_mfe].reset_index(drop=True)
    # print(candidates_sorted.head(4))

    if initial:
        no_of_guides = no_of_grnas_first
    else:
        no_of_guides = no_of_grnas_subsequent
    # print('sleeping')
    # time.sleep(5)
    output_duplexes = []
    for i in range(len(mfe_filtered['Guide_name'])):
        if len(output_duplexes) == no_of_guides:
            break
        if any(mfe_filtered['Guide_name'].loc[i] in l for l in output_duplexes):
            print(f'Guide {i} has a preferred alternative binding site.')
            pass
        else:
            candidate = mfe_filtered[['Guide_name', 'mDock']].loc[i].to_list()
            gIndex = get_index(messenger_sequence, guides_dict[candidate[0]].seq, candidate[1])
            if gIndex >= min_anchor:
                candidate.append(messenger)
                candidate.append(gIndex)
                output_duplexes.append(candidate)

    return output_duplexes
