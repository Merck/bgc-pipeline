#!/usr/bin/env python
# David Prihoda
# Get CSV of candidates that are unique based on their sequence of Pfam IDs (unique 'candidate_hash' column).

import argparse
import pandas as pd

def aggregate_unique(group):
    """
    Create a single unique candidate from its different duplications found in different locations and contigs.

    :param group: DataFrame with duplications of a single candidate, all with the same candidate_hash and sequence of pfam_ids
    :return: single unique candidate represented by a Series
    """
    contigs = group['contig_id'].unique()
    ids = group['candidate_id'].unique()
    return pd.Series({
        'num_contigs': len(contigs),
        'num_occurences': len(ids),
        'avg_prediction': group['avg_prediction'].mean(),
        'avg_nucl_length': group['nucl_length'].mean(),
        'avg_num_proteins': group['num_proteins'].mean(),
        'num_bio_domains': group['num_bio_domains'].iloc[0],
        'num_all_domains': group['num_all_domains'].iloc[0],
        'pfam_ids': group['pfam_ids'].iloc[0]
    })

def get_unique_candidates(cands):
    """
    Get DataFrame of candidates that are unique based on their sequence of Pfam IDs (unique 'candidate_hash' column).

    :param cands: DataFrame of non-unique candidates.
    :return: DataFrame of unique candidates and their number of occurences and unique contigs.
    """
    return cands.groupby('candidate_hash').apply(aggregate_unique)[[
        'num_contigs', 'num_occurences', 'avg_prediction', 'avg_nucl_length', 'avg_num_proteins', 'num_bio_domains', 'num_all_domains', 'pfam_ids'
    ]].sort_values(by='num_occurences', ascending=False)

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Target model candidate csv file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    cands = pd.read_csv(options.input)

    unique_cands = get_unique_candidates(cands)

    unique_cands.to_csv(options.output, index=True)
    print('Saved {} unique candidates to: {}'.format(len(unique_cands), options.output))

