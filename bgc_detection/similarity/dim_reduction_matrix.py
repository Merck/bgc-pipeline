#!/usr/bin/env python
# David Prihoda
# Run t-Distributed Stochastic Neighbor Embedding (t-SNE) dimensionality reduction


import argparse
import pandas as pd
from sklearn import manifold
import itertools

def create_dim_reduction_model(num_components=2, tsne_perplexity=30.0, init='random', metric='euclidean'):
    """
    Create model for t-Distributed Stochastic Neighbor Embedding (t-SNE) dimensionality reduction.
    Multicore TSNE will be used for euclidean metric and random initialization - should run faster and in parallel on LINUX systems.

    :param num_components: Reduced output dimensionality
    :param tsne_perplexity: t-SNE perplexity, influences layout (see sklearn manifold.TSNE)
    :param init: How to initialize the embedding (see sklearn manifold.TSNE)
    :param metric: Distance metric to use
    :return: Dimensionality reduction model
    """
    model_class = manifold.TSNE
    model_args = {
        'n_components': num_components,
        'init': init,
        'perplexity': tsne_perplexity,
        'metric': metric,
        'random_state': 0
    }
    if metric == 'euclidean' and init == 'random':
        try:
            from MulticoreTSNE import MulticoreTSNE
            import multiprocessing
            model_class = MulticoreTSNE
            model_args['n_jobs'] = multiprocessing.cpu_count()
            print('Using multicore t-SNE with {} cores.'.format(model_args['n_jobs']))
        except:
            print('MulticoreTSNE module not present, falling back to sklearn implementation.')
    return model_class(**model_args)

def dim_reduction_matrix(features, model):
    """
    Reduce dimensionality of given feature matrix
    :param features: DataFrame/matrix with numeric columns
    :param model: Model to use for reduction
    :return: Reduced DataFrame with columns starting from 'a', 'b', etc...
    """
    Y = model.fit_transform(features)
    return pd.DataFrame(Y, index=features.index, columns=list('abcdefghijklmnopqrstuvwxyz'[:Y.shape[1]]))

def convert_feature_metric(metric, features):
    if metric == 'similarity':
        metric = 'precomputed'
        features = (1 / features).clip(upper=1e6)
        print(features.describe())
    return metric, features

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True, action='append',
                      help="Feature matrix CSV/pickle file.", metavar="FILE")

    parser.add_argument("--index", dest="input_index", action='append', default=[],
                      help="Column in CSV feature file to use as index.", metavar="FILE")

    parser.add_argument("--num-components", dest="num_components", default=2, type=int,
                      help="Number of final components (dimensions).", metavar="INT")

    parser.add_argument("--metric", dest="metric", default='euclidean',
                      help="Metric used to compute similarity.", metavar="STRING")

    parser.add_argument("--tsne-perplexity", dest="tsne_perplexity", default=30.0, type=float,
                      help="Perplexity setting of t-SNE.", metavar="FLOAT")

    parser.add_argument("--limit", dest="limit", default=None, type=int,
                      help="Sample a limited number of vectors.", metavar="INT")

    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    all_features = []
    for input_path, input_index in itertools.zip_longest(options.input, options.input_index):
        if input_path.endswith('.csv'):
            if not input_index:
                input_index = 'contig_id'
            features = pd.read_csv(input_path)
            features.set_index(input_index, inplace=True)
        else:
            features = pd.read_pickle(input_path)
        if options.limit and len(features) > options.limit:
            features = features.sample(options.limit)
        print('Got {} feature vectors from {}'.format(features.shape, input_path))
        all_features.append(features)

    all_features: pd.DataFrame = pd.concat(all_features).fillna(0)

    metric, all_features = convert_feature_metric(options.metric, all_features)

    model = create_dim_reduction_model(
        num_components=options.num_components,
        tsne_perplexity=options.tsne_perplexity,
        metric=metric
    )
    print('Reducing feature space from {} to {}...'.format(all_features.shape, (all_features.shape[0], options.num_components)))
    matrix = dim_reduction_matrix(all_features, model)
    matrix.to_csv(options.output, index=True)
    print('Saved {} matrix to: {}'.format(matrix.shape, options.output))

