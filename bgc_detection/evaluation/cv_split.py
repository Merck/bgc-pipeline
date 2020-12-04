#!/usr/bin/env python
# David Prihoda
# Split samples from multiple Domain CSV files using n-fold cross validation
# Produces a train and test Domain CSV file for each split, along with a splits.csv file with metadata

import merge_split_samples
import argparse
import numpy as np
import pandas as pd
import os

def cv_samples(samples, result_path, n_folds, random_seed=0):
    """
    Save n_folds splits from given list of DataFrames (each DataFrame represents one sample, a BGC or non-BGC, each sample and has multiple rows - protein domains)
    Produces a train and test Domain CSV file for each split, along with a splits.csv file with metadata.

    :param samples: List of Domain CSV DataFrames
    :param result_path: Path to output directory
    :param n_folds: Number of folds (splits) to create
    :param random_seed: Random seed used for splitting
    """
    meta = []
    print('Random seed', random_seed)

    np.random.seed(random_seed)
    merged_splits = merge_split_samples.merged_split(
        samples,
        n_folds,
        shuffle_train=False,
        shuffle_test=True
    )
    merged_splits = list(merged_splits)

    i = 0
    for split_domains in merged_splits:
        i += 1
        name = 'fold{}'.format(i)
        label = 'Fold {}'.format(i)
        train_domains, test_domains = split_domains
        train_csv_path = os.path.join(result_path, name+'.train.csv')
        train_domains.to_csv(train_csv_path, index=False)
        print('Saved LCO {} train sequence to: {}'.format(name, train_csv_path))
        test_csv_path = os.path.join(result_path, name+'.test.csv')
        test_domains.to_csv(test_csv_path, index=False)
        print('Saved LCO {} test sequence to: {}'.format(name, test_csv_path))
        meta.append({
            'label': label,
            'name': name,
            'group': name
        })

    meta_csv_path = os.path.join(result_path, 'splits.csv')
    meta = pd.DataFrame(meta)
    meta.to_csv(meta_csv_path, index=False)
    print('Saved splits meta file to: {}'.format(meta_csv_path))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--maxevalue", dest="maxevalue", required=True, type=float,
                        help="Maximum domain independent e-value.", metavar="FLOAT")
    parser.add_argument("--random-seed", dest="random_seed", required=False, type=int, default=0,
                        help="Random seed used to shuffle the samples.", metavar="INT")
    parser.add_argument("-f", "--folds", dest="folds", required=True, type=int,
                        help="Number of folds.", metavar="INT")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output samples folder path.", metavar="FILE")
    parser.add_argument("--randomize-label-column", dest="randomize_label_column",
                        help="Label column to replace with random value (used for CV baseline testing).", metavar='STR')
    parser.add_argument(dest='samples', nargs='+',
                        help="Paths to CSV domain files to evaluate on.", metavar="SAMPLES")
    options = parser.parse_args()

    all_samples = []
    for path in options.samples:
        domains = pd.read_csv(path)
        domains = domains[domains['evalue'] < options.maxevalue]
        samples = [s for i, s in domains.groupby('contig_id')]
        print('Loaded {} samples and {} domains from {}'.format(len(samples), len(domains), path))
        if options.randomize_label_column:
            for sample in samples:
                sample[options.randomize_label_column] = np.random.randint(0, 2)
        all_samples += samples

    print('Output will be saved to {}/'.format(os.path.abspath(options.output)))
    os.makedirs(options.output)

    cv_samples(
        samples=all_samples,
        result_path=options.output,
        n_folds=options.folds,
        random_seed=options.random_seed
    )
