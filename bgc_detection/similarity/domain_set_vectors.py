#!/usr/bin/env python
# David Prihoda
# Create domain set (one-hot-encoding vector) representation from a Domain CSV file
# Will produce a single one-hot-encoding vector for each sample sequence (unique contig_id)

import argparse
import pandas as pd
import numpy as np

def domain_set_vectors(domains, limit=None, ratios=False):
    """
    Create one-hot-encoding representation of protein domain sequences. Will produce a column for each unique pfam ID where 0 = not present, 1 = present.
    :param domains: Domain DataFrame with contig_id and pfam_id column
    :param limit: Only add given number of most common pfam columns. Use None for no limit.
    :param ratios: Return ratio instead of binary values (number of occurences of given pfam divided by number of total pfam domains in sample)
    :return: one-hot-encoding representation of protein domain sequences, each protein domain sequence is represented by a single vector
    """
    print('Computing pfam counts from {} domains...'.format(len(domains)))
    counts = domains[['contig_id', 'pfam_id']].pivot_table(index='contig_id', columns='pfam_id',
                                                           aggfunc=np.count_nonzero).fillna(0)

    if ratios:
        counts = counts.divide(counts.sum(axis=1).values, axis=0)
    else:
        counts = counts.astype(bool).astype(int)
    if limit:
        totals = counts.sum().sort_values(ascending=False).iloc[:limit]
        print('Selecting top {} most frequent pfams:'.format(limit))
        print(totals)
        counts = counts[totals.index]
    print('Converting to sparse format...')
    return counts.to_sparse(fill_value=0)


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Domain csv file.", metavar="FILE")
    parser.add_argument("-l", "--limit", dest="limit", required=False, default=None, type=int,
                        help="Only select top N most frequent pfams.", metavar="INT")
    parser.add_argument("-r", "--ratios", dest="ratios", action='store_true', default=False,
                        help="Produce ratios instead of binary values.")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    domains = pd.read_csv(options.input)

    vectors = domain_set_vectors(domains, limit=options.limit, ratios=options.ratios)

    print('Saving pickle of {} domain set vectors to: {}'.format(len(vectors), options.output))
    vectors.to_pickle(options.output)
    print('Done.')

