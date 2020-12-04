#!/usr/bin/env python
# David Prihoda
# Plot confusion matrix from a given Domain CSV prediction file
# and prediction threshold defined by the TPR or FPR values to be achieved

import argparse
import pandas as pd
import numpy as np
from sklearn import metrics
import seaborn as sns
import matplotlib.pyplot as plt

def format_confusion_matrix(matrix, classes, figsize=(8, 4.5), title='', **kwargs):
    """
    Plot confusion matrix using provided sklearn metrics.confusion_matrix values
    :param matrix: Numpy array with values from sklearn metrics.confusion_matrix
    :param classes: Tuple for (negative, positive) class labels
    :param figsize: Figure size
    :param title: Figure title
    :param kwargs: Arguments to pass to the plot function
    :return: figure with confusion matrix plot
    """
    tn, fp, fn, tp = matrix.ravel()
    labels = np.array([
        [
            'TP = {}'.format(tp),
            'FP = {}'.format(fp),
            'Precision = {:.2f}%'.format(tp / (tp + fp) * 100)
        ],
        [
            'FN = {}'.format(fn),
            'TN = {}'.format(tn),
            ''
        ],
        [
            'TPR = {:.2f}%'.format(tp / (tp + fn) * 100),
            'FPR = {:.2f}%'.format(fp / (fp + tn) * 100),
            ''
        ]
    ])

    #sns.set_style("ticks", {"xtick.major.size": 0, "ytick.major.size": 0})
    #sns.set(font_scale=1.2)

    columns = ['Labelled ' + classes[1], 'Labelled ' + classes[0], '']
    index = ['Predicted ' + classes[1], 'Predicted ' + classes[0], '']
    vals = np.array([[tp, fp, 0], [fn, tn, 0], [0, 0, 0]])
    template = pd.DataFrame(vals, index=index, columns=columns)
    print(template)

    vmax = np.sum(vals)
    cmap = sns.blend_palette(['white', '#0066cc'], as_cmap=True)
    fig, ax = plt.subplots(1,1,figsize=figsize)
    sns.heatmap(template, ax=ax, annot=labels, fmt='', vmax=vmax, cbar=False, cmap=cmap, linewidths=1, **kwargs)
    ax.xaxis.tick_top()
    plt.suptitle(title, fontsize=13)
    plt.yticks(rotation=0)
    ax.tick_params(labelsize=15)
    fig.tight_layout()
    fig.subplots_adjust(top=0.77)

    return fig

def get_threshold(true_values, predictions, target_fpr=None, target_tpr=None):
    """
    Calculate threshold that should be used a given FPR or TPR value, based on given true values and predictions.
    Can be seen as a horizontal or vertical cut of a ROC curve
    :param true_values: Series of true values
    :param predictions: Series of predictions
    :param target_fpr: Target TPR to be achieved (or None to ignore)
    :param target_tpr: Target FPR to be achieved (or None to ignore)
    :return: threshold that should be used a given FPR or TPR value
    """
    if target_fpr is None and target_tpr is None:
        raise AttributeError('Specify one of TPR and FPR')
    if target_fpr and target_tpr:
        raise AttributeError('Specify only one of TPR and FPR')
    prev_threshold = None
    fprs, tprs, thresholds = metrics.roc_curve(true_values, predictions)
    for fpr, tpr, threshold in zip(fprs, tprs, thresholds):
        if target_fpr is not None and fpr > target_fpr:
            break
        if target_tpr is not None and tpr > target_tpr:
            break
        prev_threshold = threshold
    if not prev_threshold:
        raise AttributeError('Target FPR or TPR not achievable')

    return prev_threshold

def confusion_matrix(true_values, predictions, threshold, title, **kwargs):
    """
    Plot confusion matrix from a given Domain CSV prediction file and prediction threshold
    :param true_values: Series of true values
    :param predictions: Series of predictions
    :param threshold: Inclusive prediction threshold to use
    :param title: Plot title
    :param kwargs: Additional arguments for plot function
    :return: Figure with confusion matrix
    """
    matrix = metrics.confusion_matrix(true_values, predictions >= threshold)

    title = title + ' (threshold {:.5f})'.format(threshold)
    print(title)
    return format_confusion_matrix(matrix, title=title, classes=['non-BGC', 'BGC'], **kwargs)

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("--fpr", dest="fpr", default=None, type=float,
                      help="Target validation FPR used to select threshold.", metavar="FLOAT")
    parser.add_argument("--tpr", dest="tpr", default=None, type=float,
                      help="Target validation TPR used to select threshold.", metavar="FLOAT")
    parser.add_argument("--threshold", dest="threshold", required=False, type=float,
                      help="Prediction threshold to select domains.", metavar="INT")
    parser.add_argument("-t", "--title", dest="title", required=True,
                      help="Confusion matrix plot title.", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")
    parser.add_argument(dest='predictions', nargs='+',
                        help="Paths to CSV prediction files.", metavar="FILE")

    options = parser.parse_args()

    predictions = [pd.read_csv(path) for path in options.predictions]
    merged_true_values = np.concatenate([p['in_cluster'] for p in predictions])
    merged_predictions = np.concatenate([p['prediction'] for p in predictions])

    if options.threshold:
        threshold = options.threshold
    elif options.fpr or options.tpr:
        threshold = get_threshold(merged_true_values, merged_predictions, target_fpr=options.fpr,
                                  target_tpr=options.tpr)
    else:
        raise AttributeError('Specify either threshold or target TPR/FPR')

    fig = confusion_matrix(merged_true_values, merged_predictions, threshold, title=options.title)

    fig.savefig(options.output, dpi=100)
    print('Saved plot to {}'.format(options.output))