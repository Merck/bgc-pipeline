#!/usr/bin/env python
# David Prihoda
# Evaluate prediction results of multiple test sequences, produce per-group and averaged figures
# Used for cross-validation and leave-class-out validation

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from evaluation import evaluation_plots
from pipeline import PipelineWrapper
import os
import matplotlib

def evaluate_splits(model_folders, splits, title, result_path, figsize=(5, 5)):
    """

    :param model_folders: Paths to model folders, each folder should contain a config.json file with model config
    and splitN.test.csv and splitN.train.csv for each split in splits.csv
    :param splits: DataFrame with paths to train and test split files
    :param title: Plot title
    :param result_path: Plot output folder path
    :param figsize: Size of all figures
    :return:
    """
    os.makedirs(result_path, exist_ok=True)
    splits_by_group = splits.groupby('group')
    num_splits = len(splits_by_group)
    num_models = len(model_folders)
    num_output_rows = num_models * num_splits
    num_roc_rows = 2
    num_roc_columns = int(np.ceil(num_splits / num_roc_rows))

    probability_fig, probability_ax = plt.subplots(num_output_rows, 1, figsize=(100, num_output_rows * 1.5))

    samples_fig, samples_ax = plt.subplots(num_roc_rows, num_roc_columns, figsize=(num_roc_columns * figsize[0], num_roc_rows * figsize[1]))
    if num_roc_columns == 1:
        samples_ax = [samples_ax]

    mean_fig, mean_ax = plt.subplots(1, 1, figsize=figsize)

    for p, path in enumerate(model_folders):
        pipeline = PipelineWrapper.from_config(os.path.join(path, 'config.json'), meta_only=True)
        print('='*80)
        print(pipeline.label)
        print('='*80)
        # Train and validate each model
        predictions = []
        true_outputs = []
        split_no = 0
        for group_name, group_splits in splits_by_group:
            group_label = group_splits['label'].iloc[0]
            test_domains: pd.DataFrame = pd.concat([pd.read_csv(os.path.join(path, split_name+'.test.csv')) for split_name in group_splits['name']])
            prediction = test_domains['prediction']
            predictions.append(prediction)

            true_output = test_domains['in_cluster']
            true_outputs.append(true_output)

            # Plot sample ROC
            print(group_name)
            evaluation_plots.plot_roc_curve(
                true_output,
                prediction,
                ax=samples_ax[split_no//num_roc_columns][split_no % num_roc_columns],
                title=group_label,
                label=pipeline.label,
                color=pipeline.color
            )

            # Plot sample output
            ax = probability_ax[split_no*num_models + p]
            prob_title = group_name + ': ' + pipeline.label
            df = test_domains.reset_index(drop=True)
            df['in_cluster'].plot(kind='area', ax=ax, color='black', lw=0.5, alpha=0.2, label=None)
            df['in_cluster'].plot(ax=ax, color='black', lw=0.5, alpha=0.9, label='true_output')
            df['prediction'].plot(ax=ax, title=prob_title, lw=0.5, alpha=0.7, label=pipeline.label, color=pipeline.color)
            ax.set_ylim([-0.05, 1.05])

            split_no += 1

        print('-'*80)
        print('Mean ROC:')

        evaluation_plots.plot_roc_curve(
            np.concatenate(true_outputs),
            np.concatenate(predictions),
            ax=mean_ax,
            title=title,
            lw=1,
            label=pipeline.label,
            color=pipeline.color
        )

        print('-' * 80)

        mean_fig_path = os.path.join(result_path, 'roc_mean.png')
        mean_fig.savefig(mean_fig_path, dpi=150)
        print('Saved mean ROC plot to: ', mean_fig_path)
        mean_fig_path_pdf = os.path.join(result_path, 'roc_mean.pdf')
        mean_fig.savefig(mean_fig_path_pdf, bbox_inches='tight')
        print('Saved mean ROC plot to: ', mean_fig_path_pdf)

        samples_fig_path = os.path.join(result_path, 'roc_samples.png')
        samples_fig.tight_layout(h_pad=1.8)
        samples_fig.savefig(samples_fig_path)
        print('Saved per-sample ROC plot to: ', samples_fig_path)
        samples_fig_path_pdf = os.path.join(result_path, 'roc_samples.pdf')
        samples_fig.savefig(samples_fig_path_pdf, bbox_inches='tight')
        print('Saved per-sample ROC plot to: ', samples_fig_path_pdf)

        probability_fig_path = os.path.join(result_path, 'probability.png')
        probability_fig.tight_layout()
        probability_fig.savefig(probability_fig_path)
        print('Saved per-sample predictions plot to: ', probability_fig_path)


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Splits meta CSV file path.", metavar="FILE")
    parser.add_argument("-t", "--title", dest="title", required=True,
                        help="Plot title.", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output folder path.", metavar="FILE")
    parser.add_argument("--size", dest="size", required=False, default=5, type=float,
                      help="Figure size.", metavar="FLOAT")
    parser.add_argument(dest='models', nargs='+',
                        help="Paths to model folders with configs and predictions.", metavar="SAMPLES")

    options = parser.parse_args()

    meta = pd.read_csv(options.input, low_memory=False)

    splits_folder = os.path.dirname(options.input)

    font = {'family': 'Arial', 'size': 12}
    matplotlib.rc('font', **font)
    matplotlib.rc('legend', fontsize=11, handlelength=2)

    evaluate_splits(
        model_folders=options.models,
        splits=meta,
        title=options.title,
        result_path=options.output,
        figsize=(options.size, options.size)
    )



