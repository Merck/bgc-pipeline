#!/usr/bin/env python
# David Prihoda
# Create pfam2vec mean representation from a Domain CSV file
# Will produce a single pfam2vec-mean vector for each sample sequence (unique contig_id)

import argparse
import pandas as pd
import word2vec

def pfam2vec_mean_vectors(domains, model, quantiles=None):
    """
    Create pfam2vec-mean representation of protein domain sequences.
    :param domains: Domain DataFrame with contig_id and pfam_id column
    :param model: Trained word2vec (pfam2vec) model.
    :param quantiles: Use given quantile values instead of mean
    :return: pfam2vec-mean representation of protein domain sequences, each protein domain sequence is represented by a single vector - mean of all pfam2vec vectors
    """

    pfam_vectors = pd.DataFrame(model.vectors)
    pfam_vectors['pfam_id'] = model.vocab
    pos_features = domains[['pfam_id','contig_id']].merge(pfam_vectors, on='pfam_id', how='inner')
    if quantiles:
        print('Using quantiles: {}'.format(quantiles))
        quantile_vectors = pos_features.groupby('contig_id').quantile(quantiles).unstack()
        quantile_vectors.columns = ['{}q{}'.format(dim, q) for dim, q in quantile_vectors.columns.values]
        return quantile_vectors
    else:
        print('Using mean value.')
        return pos_features.groupby('contig_id').mean()

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Domain csv file.", metavar="FILE")
    parser.add_argument("-p", "--pfam2vec", dest="pfam2vec", required=True,
                      help="Pfam2vec model file.", metavar="FILE")
    parser.add_argument("-q", "--quantile", dest="quantiles", required=False, action='append', type=float,
                      help="Use quantile(s). Repeat argument for more quantiles.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    domains = pd.read_csv(options.input)
    model = word2vec.load(options.pfam2vec)

    vectors = pfam2vec_mean_vectors(domains, model, options.quantiles)

    print('Saving {} pfam2vec vectors to: {}'.format(len(vectors), options.output))
    vectors.to_pickle(options.output)
    print('Done.')
