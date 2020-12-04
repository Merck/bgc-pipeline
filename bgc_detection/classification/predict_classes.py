#!/usr/bin/env python
# David Prihoda
# Predict classes using a vector representation of samples and a trained multilabel classifier

import pandas as pd
import numpy as np
import pickle
import argparse

def predict_classes(samples, model, add_classes_list=True):
    # Set missing columns to 0
    if not hasattr(model, 'input_columns'):
        raise AttributeError('Trained model does not contain the "input_columns" attribute.')
    if not hasattr(model, 'label_columns'):
        raise AttributeError('Trained model does not contain the "label_columns" attribute.')

    missing_columns = set(model.input_columns).difference(samples.columns)
    for col in missing_columns:
        samples[col] = 0
    print('Missing columns:\n{}'.format(sorted(list(missing_columns))))
    print('Warning: Setting {} missing columns to 0'.format(len(missing_columns)))
    samples = samples[model.input_columns]

    results = np.array([r[:,1] for r in model.predict_proba(samples.values)]).transpose()
    predictions = pd.DataFrame(results, index=samples.index, columns=model.label_columns)
    if add_classes_list:
        predictions['classes'] = [';'.join(model.label_columns[x >= 0.5]) for x in results]

    return predictions


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Sample vectors pickle file path (indexed by sample ID).", metavar="FILE")
    parser.add_argument("-m", "--model", dest="model", required=True,
                      help="Trained model pickle file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output predictions CSV file path.", metavar="FILE")

    options = parser.parse_args()

    samples = pd.read_pickle(options.input)

    with open(options.model, 'rb') as f:
        model = pickle.load(f)

    predictions = predict_classes(samples, model)
    predictions.to_csv(options.output)
    print('Saved {} predictions to: {}'.format(len(predictions), options.output))

