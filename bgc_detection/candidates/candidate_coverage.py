#!/usr/bin/env python
# David Prihoda
# Calculate coverage of BGCs by a DataFrame of BGC Candidates

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def get_single_contig_coverage(a_cands, b_cands):
    """
    Get coverage of each BGC candidate in a_cands by BGC candidates in b_cands,
    where all candidates come from the same contig.
    :param a_cands: Reference DataFrame of BGC candidates (from a single contig)
    :param b_cands: Compared DataFrame of BGC candidates (from a single contig)
    :return: row of each BGC candidate in a_cands with 'coverage' column that defines
            fractional coverage by overlapping BGC candidates in b_cands
    """
    if b_cands is None:
        remaining_cands = []
    else:
        remaining_cands = list(b_cands.reset_index(drop=True).iterrows())
    # Create binary mask based on longest canidate length
    max_len = int((a_cands['nucl_end'] - a_cands['nucl_start'] + 1).max())
    mask = np.zeros(max_len)
    # For each A candidate
    coverages = []
    for c, cand in a_cands.iterrows():
        # For each suitable candidate from other model
        cand_start = int(cand['nucl_start']) - 1
        cand_end = int(cand['nucl_end'])
        cand_len = cand_end - cand_start
        #print('Cand {}: {}-{} (len {})'.format(c, cand_start, cand_end, cand_len))
        any_exact = False
        max_covered = 0
        for i, other in remaining_cands:
            other_start = int(other['nucl_start']) - 1
            other_end = int(other['nucl_end'])
            other_len = other_end - other_start
            # No overlap anymore
            if other_start > cand_end:
                continue
            # No overlap yet
            if other_end < cand_start:
                # Discard all previous candidates up to current one
                continue
            # Exact match
            if other_start == cand_start and other_end == cand_end:
                any_exact = True
            # Start and end coordinates relative from cand_start
            overlap_start = max(other_start, cand_start) - cand_start
            overlap_end = min(other_end, cand_end) - cand_start
            overlap_length = overlap_end - overlap_start
            mask[overlap_start:overlap_end] = 1
            max_covered = max(max_covered, overlap_length / other_len)

        num_covered = sum(mask[:cand_len])
        mask[:cand_len] = 0

        #print('overlap {}/{} = {}'.format(num_covered, cand_len, num_covered / cand_len))
        coverage = pd.Series(
            [num_covered / cand_len, any_exact, max_covered],
            ['coverage', 'any_exact', 'max_covered']
        ).append(cand)

        if 'model' in coverage:
            del coverage['model']
        coverages.append(coverage)
    return coverages

def get_coverage(a_cands, b_cands):
    """
    Get coverage of each BGC candidate in a_cands by BGC candidates in b_cands,
    where each candidate can be found in a different contig as indicated by the 'contig_id' column.
    :param a_cands: Reference DataFrame of BGC candidates
    :param b_cands: Compared DataFrame of BGC candidates
    :return: row of each BGC candidate in a_cands with 'coverage' column that defines
            fractional coverage by overlapping BGC candidates in b_cands
    """
    a_grouped = a_cands.groupby('contig_id')
    b_grouped = b_cands.groupby('contig_id')

    coverages = []
    # Get coverage separately for all contigs
    for contig_id in a_grouped.groups:
        #print(contig_id)
        a_contig_cands = a_grouped.get_group(contig_id)
        b_contig_cands = b_grouped.get_group(contig_id) if contig_id in b_grouped.groups else None
        coverages += get_single_contig_coverage(a_contig_cands, b_contig_cands)
    return pd.DataFrame(coverages)


def plot_coverage_hist(coverage, title, label, **kwargs):
    """
    Plot histogram of coverage by model
    :param coverage: DataFrame with BGC candidates and their 'coverage' column and 'model' column
    :param title: Plot title
    :param label: Plot x-axis label
    :param kwargs: Arguments to pass to the histogram plot function
    """
    cols = len(coverage['model'].unique())
    axes = coverage[['coverage', 'model']].hist(by='model', bins=25, figsize=(cols * 3, 2.7), layout=(1, cols),
                                                sharey=True, **kwargs)
    axes[0].set_ylabel('# BGCs')
    plt.suptitle(title)
    plt.tight_layout()
    plt.subplots_adjust(top=0.77)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_xlabel(label)
        ax.set_xticklabels(['{:.0f}%'.format(x * 100) for x in ax.get_xticks()])


def plot_coverage_boxplot(coverage, title, label, **kwargs):
    """
    Plot boxplot of coverage by model
    :param coverage: DataFrame with BGC candidates and their 'coverage' column and 'model' column
    :param title: Plot title
    :param label: Plot x-axis label
    :param kwargs: Arguments to pass to the boxplot function
    """
    cols = len(coverage['model'].unique())
    ax = coverage[['coverage', 'model']].boxplot(by='model', figsize=(cols * 0.7+1, 2.7), **kwargs)
    plt.suptitle(title)
    plt.tight_layout()
    plt.xticks(rotation=90)
    plt.subplots_adjust(top=0.80)
    ax.set_ylabel(label)
    ax.set_yticklabels(['{:.0f}%'.format(x * 100) for x in ax.get_yticks()])


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Target model candidate csv file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")
    parser.add_argument(dest='candidates', nargs='+',
                        help="Paths to other models' candidate files.", metavar="FILE")

    options = parser.parse_args()

    target_cands = pd.read_csv(options.input)
    other_cands: pd.DataFrame = pd.concat([pd.read_csv(path) for path in options.candidates])

    coverage = get_coverage(target_cands, other_cands)

    coverage.to_csv(options.output, index=False)
    print('Saved {} candidates to: {}'.format(len(coverage), options.output))

