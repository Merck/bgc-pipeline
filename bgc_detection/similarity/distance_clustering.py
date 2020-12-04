#!/usr/bin/env python
# David Prihoda
# Cluster a numeric DataFrame using different clustering methods

import argparse
import pandas as pd
import itertools
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
import numpy as np

def distance_clustering(features, n_clusters, method='kmeans', metric='euclidean'):
    """
    Cluster a DataFrame with numeric columns into n_cluster clusters using given algorithm
    :param features: Feature matrix (DataFrame)
    :param n_clusters: Number of clusters
    :param method: Clustering method
    :param metric: Distance metric
    :return: Cluster labels, one for each row in input DataFrame.
    """
    if method == 'kmeans':
        if metric != 'euclidean':
            raise ValueError('Only euclidean metric is allowed for KMeans.')
        model = KMeans(n_clusters=n_clusters, random_state=0)
    elif method == 'agglomerative':
        model = AgglomerativeClustering(n_clusters=n_clusters, affinity=metric, linkage='ward')
    elif method == 'dbscan':
        model = DBSCAN(metric=metric)
    else:
        raise ValueError('Invalid clustering method {}.'.format(method))
    predictions = model.fit_predict(features)
    n_clusters = np.max(predictions) + 1
    cluster_labels = ['cluster{}'.format(i+1) for i in range(0, n_clusters)]
    sample_labels_list = [cluster_labels[idx] for idx in predictions]
    sample_labels = pd.Series(sample_labels_list, index=features.index, name='cluster')
    return sample_labels


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True, action='append',
                      help="Feature matrix CSV/pickle file.", metavar="FILE")

    parser.add_argument("--index", dest="input_index", action='append', default=[],
                      help="Column in feature file to use as index (CSV only).", metavar="FILE")

    parser.add_argument("--num-clusters", dest="num_clusters", required=False, type=int,
                      help="Number of clusters.", metavar="INT")

    parser.add_argument("--metric", dest="metric", default='euclidean',
                      help="Metric used to compute similarity.", metavar="STRING")

    parser.add_argument("--method", dest="method", default='kmeans',
                      help="Clustering method (kmeans, agglomerative).", metavar="STRING")

    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    all_features = []
    for input_path, input_index in itertools.zip_longest(options.input, options.input_index):
        if input_path.endswith('.csv'):
            if not input_index:
                input_index = 'contig_id'
            features = pd.read_csv(input_path).set_index(input_index)
        else:
            features = pd.read_pickle(input_path)
        print('Got {} feature vectors from {}'.format(features.shape, input_path))
        all_features.append(features)

    all_features: pd.DataFrame = pd.concat(all_features).fillna(0)

    print('Clustering on {} features:'.format(len(all_features.columns)))
    print(all_features.describe())
    print('Clustering...')
    clusters = distance_clustering(all_features, n_clusters=options.num_clusters, method=options.method, metric=options.metric)

    print('-' * 80)
    print(clusters.head())

    print('-'*80)
    print(clusters.value_counts())

    clusters.to_csv(options.output, index=True, header=True)
    print('Saved {} cluster labels to: {}'.format(len(clusters), options.output))

