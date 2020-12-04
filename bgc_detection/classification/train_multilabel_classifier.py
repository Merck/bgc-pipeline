#!/usr/bin/env python
# David Prihoda
# Train multilabel classifier (zero or more classes)
# using vector representation of samples and binary representation of labels

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import pickle
import argparse

def train_multilabel_classifier(samples, labels):
    """
    Train model on given vector representation of samples and binary representation of labels
    Only samples found in both datasets are used.

    :param samples: DataFrame of vectors, one for each sample, indexed by sample ID
    :param labels: DataFrame of binary label vectors, one for each sample, indexed by sample ID
    :return: trained classifier
    """
    # Get rows found in both dataframes and align the axis
    labels = labels.loc[samples.index].dropna()
    samples = samples.loc[labels.index]
    print('Got {} labelled samples'.format(len(samples)))

    print('Training model with {} classes: {}'.format(len(labels.columns.values), labels.columns.values))
    model = RandomForestClassifier(n_estimators=100, random_state=0)
    model.fit(samples.values, labels.values)
    model.input_columns = samples.columns.values
    model.label_columns = labels.columns.values

    return model


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Training sample vectors pickle file path (indexed by sample ID).", metavar="FILE")
    parser.add_argument("-c", "--classes", dest="classes", required=True,
                      help="Sample labels CSV file path (with contig_id column corresponding to sample IDs).", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output trained model file path.", metavar="FILE")

    options = parser.parse_args()

    print('Reading samples and labels...')
    samples = pd.read_pickle(options.input)
    labels = pd.read_csv(options.classes).set_index('contig_id')

    model = train_multilabel_classifier(samples, labels)

    with open(options.output, 'wb') as f:
        pickle.dump(model, f)
    print('Saved model to: {}'.format(options.output))

