#!/usr/bin/env python
# David Prihoda
# Plot boxplots of different statistics for given BGC candidate CSV files

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

def candidate_stats(cands, title=None):
    """
    Plot boxplots of different statistics for given BGC candidate DataFrame
    :param cands: DataFrame of candidates with 'model' field
    :param title: Plot title
    :return: Figure of boxplots of different statistics for given BGC candidate DataFrame
    """
    cands_by_model = cands.groupby('model')

    fig, axes = plt.subplots(2, 3, figsize=(17, 2 + 0.4 * len(cands_by_model)))
    ax = sns.barplot(0, 'model', data=cands_by_model.size().reset_index(), ax=axes[0][0], color='grey')
    ax.set_xlabel('# BGC candidates total')

    print(cands_by_model['candidate_hash'].nunique().reset_index().head())

    ax = sns.barplot('candidate_hash', 'model', data=cands_by_model['candidate_hash'].nunique().reset_index(), ax=axes[0][1], color='grey')
    ax.set_xlabel('# Unique BGC candidates')

    ax = sns.boxplot('num_proteins', 'model', ax=axes[0][2], data=cands, color='white', linewidth=1, showfliers=False)
    ax.set_xlabel('# proteins per BGC candidate')

    ax = sns.boxplot('num_all_domains', 'model', ax=axes[1][0], data=cands, color='white', linewidth=1, showfliers=False)
    ax.set_xlabel('# protein domains per BGC candidate')

    ax = sns.boxplot('num_bio_domains', 'model', ax=axes[1][1], data=cands, color='white', linewidth=1, showfliers=False)
    ax.set_xlabel('# biosynthetic protein domains per BGC candidate')

    ax = sns.boxplot('nucl_length', 'model', ax=axes[1][2], data=cands, color='white', linewidth=1, showfliers=False)
    ax.set_xlabel('nucleotide length per BGC candidate')

    fig.tight_layout()

    if title:
        fig.suptitle(title, fontsize=15)
        fig.subplots_adjust(top=0.90)

    return fig

def read_cands(paths, labels):
    """
    Read candidate CSV paths into a concatenated DataFrame with 'model' field set to corresponding label
    :param paths: List of candidate CSV paths
    :param labels: List of 'model' labels for given path
    :return: concatenated DataFrame from given candidate CSV paths with 'model' field set to corresponding label
    """
    cands = []
    for path, label in zip(paths, labels):
        c = pd.read_csv(path)
        c['model'] = label
        cands.append(c)

    cands: pd.DataFrame = pd.concat(cands)
    return cands

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")
    parser.add_argument("--candidates", dest="candidates", required=True, action='append',
                      help="Model candidates CSV file path.", metavar="FILE")
    parser.add_argument("--label", dest="labels", required=True, action='append',
                      help="Model label.", metavar="STRING")
    parser.add_argument("-t", "--title", dest="title", required=False,
                      help="Plot title.", metavar="FILE")

    options = parser.parse_args()

    cands = read_cands(options.candidates, labels=options.labels)

    fig = candidate_stats(cands, title=options.title)
    fig.savefig(options.output, dpi=200)

    print('Saved stats plot to: ', options.output)

