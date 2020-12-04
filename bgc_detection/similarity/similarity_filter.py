#!/usr/bin/env python
# David Prihoda
# Filter a Domain CSV file by comparing its feature representation to a feature representation of a different set of samples
# Will produce the original Domain CSV file where samples with distance < threshold are removed.

import argparse
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS

def get_similar_index(input_features, compare_features, metric, threshold):
    """
    Get index of input samples more distant from a second set of samples than given threshold
    :param input_features: The input samples vector representation DataFrame, one sample vector per row
    :param compare_features: The compared samples vector representation DataFrame, one sample vector per row
    :param metric: Metric to use (euclidean, cosine, ...)
    :param threshold: Return index of samples with distance is bigger than the threshold
    :return: index of input samples more distant from the second set of samples than given threshold
    """
    cols = np.intersect1d(input_features.columns, compare_features.columns)
    input_features = input_features[cols]
    compare_features = compare_features[cols]
    fn = PAIRWISE_DISTANCE_FUNCTIONS[metric]
    distances = pd.DataFrame(fn(input_features.values, compare_features.values), index=input_features.index, columns=compare_features.index)
    min_distances = distances.min(axis=1)
    print('Distance stats:')
    print(min_distances.describe())
    filtered = min_distances[min_distances >= threshold]
    print('Removed {} contigs with distance < {}'.format(len(input_features)-len(filtered), threshold))
    return filtered.index

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Input domain CSV file.", metavar="FILE")

    parser.add_argument("-f", "--feature", dest="features", required=True,
                      help="Input feature matrix pickle file.", metavar="FILE")

    parser.add_argument("-c", "--compare", dest="compare", required=True,
                      help="Feature matrix pickle file to compare against.", metavar="FILE")

    parser.add_argument("--metric", dest="metric", default='euclidean',
                      help="Metric used to compute similarity.", metavar="STRING")

    parser.add_argument("--threshold", dest="threshold", required=True, type=float,
                      help="Remove samples with distance < threshold.", metavar="FLOAT")

    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    domains = pd.read_csv(options.input).set_index('contig_id')
    input_features = pd.read_pickle(options.features)
    compare_features = pd.read_pickle(options.compare)

    print('Filtering {} contigs by similarity to {} contigs...'.format(len(input_features), len(compare_features)))
    index = get_similar_index(input_features, compare_features, metric=options.metric, threshold=options.threshold)
    orig_len = len(input_features)
    new_len = len(index)
    print('Remaining contigs: {}/{} = {:.1f}%'.format(new_len, orig_len, (new_len/orig_len)*100))
    filtered_domains = domains.loc[index]
    print('Saving {} domains to: {}'.format(len(filtered_domains), options.output))
    filtered_domains.to_csv(options.output, index=True)

