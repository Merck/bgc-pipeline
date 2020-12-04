#!/usr/bin/env python
# David Prihoda
# Calculate Levenshtein Pfam ID similarity for given Domain CSV files
# Each Pfam ID is considered a unique symbol (character)
# Will produce a symmetrical matrix

import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool
import Levenshtein


def get_sequence_similarity(first, second):
    # Sequence matcher is not symmetrical
    # return SequenceMatcher(None, first, second, autojunk=False).ratio()
    return Levenshtein.ratio(first, second)

def sequence_similarity_task(args):
    string, other_strings = args
    similarity = np.zeros(len(other_strings))
    for i, other_string in enumerate(other_strings):
        similarity[i] = get_sequence_similarity(string, other_string)
    return similarity

def get_string_from_word_list(word_list, word_index):
    return ''.join([chr(word_index[w]+200) for w in word_list])

def get_word_index(words):
    return {word: i for i, word in enumerate(set(words))}

def get_sequences_as_strings(domains, word_index):
    return domains.groupby('contig_id')['pfam_id'].apply(lambda pfam_ids: get_string_from_word_list(pfam_ids, word_index))

def sequence_similarity_matrix(a_domains, b_domains=None):
    """
    Get sequence similarity for all pairs of samples from a_domains and b_domains
    :param a_domains: Domain DataFrame with pfam_id and contig_id of samples
    :param b_domains: Domain DataFrame with pfam_id and contig_id of samples to compare with (can be equal to a_samples)
    :return: Similarity matrix indexed by samples in a_domains and columns by b_domains.
    """
    if b_domains is None:
        b_domains = a_domains
    word_index = get_word_index(np.concatenate([a_domains['pfam_id'], b_domains['pfam_id']]))
    a_strings = get_sequences_as_strings(a_domains, word_index)
    b_strings = get_sequences_as_strings(b_domains, word_index)

    print('Computing {:,} * {:,} similarities'.format(len(a_strings), len(b_strings)))

    pool = Pool()
    similarity = pool.map(sequence_similarity_task, [(string, b_strings) for string in a_strings.values])

    idx = pd.Series(a_strings.index, name='contig_id')
    return pd.DataFrame(similarity, index=idx, columns=b_strings.index)


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--compare", dest="compare", required=False,
                        help="Sequences to compare against (leave blank for self-comparison).", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file path.", metavar="FILE")
    parser.add_argument(dest='domains', nargs='+',
                        help="Domain CSV files.", metavar="FILE")

    options = parser.parse_args()

    domains: pd.DataFrame = pd.concat([pd.read_csv(path)[['contig_id', 'pfam_id']] for path in options.domains])
    to_compare = pd.read_csv(options.compare)[['contig_id', 'pfam_id']] if options.compare else None

    matrix = sequence_similarity_matrix(domains, to_compare)

    matrix.to_csv(options.output, index=True)
    print('Saved {}x{} similarity matrix to: {}'.format(matrix.shape[0], matrix.shape[1], options.output))

