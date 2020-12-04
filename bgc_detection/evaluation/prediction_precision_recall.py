#!/usr/bin/env python
# David Prihoda
# Plot Precision-Recall plot from Domain CSV with predictions

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import evaluation_plots as evaluation
import matplotlib

def prediction_precision_recall_plot(predictions_per_model, names, colors, title, figsize=(5, 5)):
    """
    Plot Precision-Recall plot from list of domain predictions (Domain DataFrames with prediction and in_cluster column)
    :param predictions_per_model: list of domain DataFrames with prediction and in_cluster column, one for each model
    :param names: List of names of models
    :param colors: List of colors of models
    :param title: Plot title
    :param figsize: Plot figure size
    :return: Plot figure.
    """
    mean_fig, mean_ax = plt.subplots(1, 1, figsize=figsize)

    for predictions, name, color in zip(predictions_per_model, names, colors):
        true_values = predictions['in_cluster']

        evaluation.plot_precision_recall_curve(
            true_values,
            predictions['prediction'],
            ax=mean_ax,
            title=title,
            lw=1,
            label=name,
            color=color
        )

    return mean_fig

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction", dest="predictions", required=True, action='append',
                      help="Prediction CSV file path.", metavar="FILE")
    parser.add_argument("--name", dest="names", required=True, action='append',
                      help="Name for given prediction.", metavar="FILE")
    parser.add_argument("--color", dest="colors", required=True, action='append',
                      help="Color for given prediction.", metavar="FILE")
    parser.add_argument("--title", dest="title", required=True,
                      help="Plot title.", metavar="FILE")
    parser.add_argument("--size", dest="size", required=False, default=5, type=int,
                      help="Figure size.", metavar="INT")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    predictions = [pd.read_csv(path) for path in options.predictions]

    font = {'family': 'Arial', 'size': 12}
    matplotlib.rc('font', **font)
    matplotlib.rc('legend', fontsize=11, handlelength=2)

    fig = prediction_precision_recall_plot(predictions, options.names, options.colors, options.title, figsize=(options.size, options.size))
    fig.savefig(options.output, dpi=150, bbox_inches='tight')

    print('Saved PR plot to: ', options.output)


