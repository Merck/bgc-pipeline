#!/usr/bin/env python
# David Prihoda
# Plot prediction probability output, producing a line plot with protein domains on x-axis and prediction value on y-axis

import argparse
import matplotlib.pyplot as plt
import pandas as pd

def count_y_clusters(y):
    prev = 0
    clusters = 0
    for val in y:
        if val == 1 and prev == 0:
            clusters += 1
        prev = val
    return clusters

def prediction_plot(predictions_per_model, names, colors):
    """
    Plot prediction probability output, producing a line plot with protein domains on x-axis and prediction and in_cluster value on y-axis
    :param predictions_per_model: list of model predictions (each model prediction is a Domain DataFrame with prediction and in_cluster column)
    :param names: list of model names
    :param colors: list of model colors
    :return: Plot figure
    """
    contig_ids = predictions_per_model[0]['contig_id'].unique()
    rows = len(contig_ids) * len(names)
    samples_fig, samples_ax = plt.subplots(rows, 1, figsize=(20, rows*2))

    i = 0
    for contig_id in contig_ids:
        for predictions, name, color in zip(predictions_per_model, names, colors):
            contig_predictions = predictions[predictions['contig_id'] == contig_id].reset_index(drop=True)
            print('Plotting', contig_id, name, len(contig_predictions))

            contig_predictions['in_cluster'].plot.area(color='black', alpha=0.1, ax=samples_ax[i])
            contig_predictions['in_cluster'].plot(color='black', lw=2, alpha=0.7, ax=samples_ax[i])
            contig_predictions['prediction'].plot(color=color, title=contig_id, label=name, ax=samples_ax[i])
            i += 1

    samples_fig.tight_layout()
    return samples_fig

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction", dest="predictions", required=True, action='append',
                      help="Prediction CSV file path.", metavar="FILE")
    parser.add_argument("--name", dest="names", required=True, action='append',
                      help="Name for given prediction.", metavar="FILE")
    parser.add_argument("--color", dest="colors", required=True, action='append',
                      help="Color for given prediction.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    predictions = [pd.read_csv(path) for path in options.predictions]

    fig = prediction_plot(predictions, options.names, options.colors)
    fig.savefig(options.output, dpi=150)

    print('Saved prediction plot to: ', options.output)


