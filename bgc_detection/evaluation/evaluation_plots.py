#!/usr/bin/env python
# David Prihoda
# Various helper functions for BGC detector evaluation

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

import numpy as np
from sklearn.metrics import precision_recall_curve, average_precision_score

def plot_roc_curve(true_values, predictions, ax=None, title='ROC', label='ROC', lw=1, add_auc=True, baseline=True, figsize=(5, 5), **kwargs):
    """
    Plot ROC curve of a single model. Can be called repeatedly with same axis to plot multiple curves.
    :param true_values: Series of true values
    :param predictions: Series of prediction values
    :param ax: Use given axis (will create new one if None)
    :param title: Plot title
    :param label: ROC curve label
    :param lw: Line width
    :param add_auc: Add AUC value to label
    :param baseline: Plot baseline that indicates performance of random model (AUC 0.5)
    :param figsize: Figure size
    :param kwargs: Additional arguments for plotting function
    :return: Figure axis
    """
    fpr, tpr, _ = roc_curve(true_values, predictions)
    roc_auc = auc(fpr, tpr)
    label_auc = label + ': {:.3f} AUC'.format(roc_auc)
    print(label_auc)
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    if baseline:
        ax.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')
    ax.plot(fpr, tpr, lw=lw, label=label_auc if add_auc else label, **kwargs)
    ax.set_title(title)
    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.legend(loc='lower right', frameon=False)
    return ax


def plot_precision_recall_curve(true_values, predictions, ax=None, title='Precision-Recall', label='PR', **kwargs):
    """
    Plot Precision-Recall curve of a single model. Can be called repeatedly with same axis to plot multiple curves.
    :param true_values: Series of true values
    :param predictions: Series of prediction values
    :param ax: Use given axis (will create new one if None)
    :param title: Plot title
    :param label: ROC curve label
    :param kwargs: Additional arguments for plotting function
    :return: Figure axis
    """
    precision, recall, _ = precision_recall_curve(true_values, predictions)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    avg_precision = average_precision_score(true_values, predictions)
    ax.step(recall, precision, where='post', label=label + ': {:.3f} AvgPrec'.format(avg_precision), **kwargs)

    ax.set_title(title)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlim([0.0, 1.0])
    ax.legend(loc='lower left', frameon=False)

    return ax


def plot_samples_roc(groups_true_values, groups_predictions, groups_titles,
                         fig = None, ax = None, columns = 5, size=4, **kwargs):
    """
    Plot multiple ROC curves for multiple groups of predictions
    :param groups_true_values: List of lists (true values)
    :param groups_predictions: List of lists (predictions)
    :param groups_titles: List of lists (plot titles)
    :param fig: Use given figure (will be created if None)
    :param ax: Use given axis (will be created if None)
    :param columns: Number of ROC axes in single row
    :param size: Size of each ROC axis
    :param kwargs: Additional plotting function arguments
    :return: Figure and axis of final plot
    """
    if fig is None:
        num_groups = len(groups_true_values)
        rows = int(np.ceil(num_groups / columns))
        fig, ax = plt.subplots(rows, columns, figsize=(columns * size, rows * size))
        # Hack: make sure ax is a two-dimensional array even if we only have one row
        if rows == 1:
            ax = [ax]

    for i, (true_values, predictions, group_title) in enumerate(zip(groups_true_values, groups_predictions, groups_titles)):
        #df = pd.DataFrame({'prediction': predictions, 'true_values': true_values}).reset_index(drop=True)
        #df.plot(figsize=(60, 2), title=group_title)
        print(group_title)
        plot_roc_curve(true_values, predictions, ax=ax[i//columns][i % columns], title=group_title, **kwargs)

    fig.tight_layout()

    return fig, ax

def set_axis_fontsize(ax, fontsize):
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_legend().get_texts() + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
