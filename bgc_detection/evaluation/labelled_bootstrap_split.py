#!/usr/bin/env python
# David Prihoda
# Create bootstrap split from Cimermanic et. al. labelled genomes
# Generates a folder with (train, test) files for each bootstrap split
# Splits the genomes based on Genome ID, where all contigs from given genome are considered together

import argparse
import numpy as np
import pandas as pd
import os
from pprint import pprint

def genome_ids_to_contig_ids(summary, genome_ids):
    contig_ids = summary[summary['Genome ID'].apply(lambda g: g in genome_ids)]['NCBI ID']
    return [c for c in contig_ids if c != '?']

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Path to domain CSV file.", metavar="FILE")
    parser.add_argument("-s", "--summary", dest="summary", required=True,
                        help="Path to labelled contig summary CSV file.", metavar="FILE")
    parser.add_argument("-n", "--number", dest="number", required=True, type=int,
                        help="Number of splits to generate.", metavar="INT")
    parser.add_argument("-r", "--train-ratio", dest="train_ratio", required=True, type=float,
                        help="Ratio of samples to use for training (rest used for testing).", metavar="FLOAT")
    parser.add_argument("--random-seed", dest="random_seed", default=1, type=int,
                        help="Random seed for splitting.", metavar="INT")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Path to output splits folder.", metavar="FILE")
    options = parser.parse_args()

    summary = pd.read_csv(options.summary, sep=';')
    domains = pd.read_csv(options.input)

    print('Output will be saved to:', options.output)
    os.makedirs(options.output)

    genome_ids = summary[summary['NCBI ID'] != '?']['Genome ID'].unique()
    print('Splitting genome IDs:', genome_ids)

    splits = []
    np.random.seed(options.random_seed)
    train_num = int(np.round(len(genome_ids) * options.train_ratio))
    for i in range(0, options.number):
        train_genomes = np.random.choice(genome_ids, train_num, replace=True)
        test_genomes = np.setdiff1d(genome_ids, train_genomes)
        train_contigs = genome_ids_to_contig_ids(summary, train_genomes)
        test_contigs = genome_ids_to_contig_ids(summary, test_genomes)

        splits.append({
            'train_genomes': ','.join(train_genomes),
            'test_genomes': ','.join(test_genomes),
            'train_contigs': ','.join(train_contigs),
            'test_contigs': ','.join(test_contigs),
            'index': i
        })
        pprint(splits[-1])

        split_train_path = os.path.join(options.output, 'split_{}_{}.csv'.format(i, 'train'))
        split_test_path = os.path.join(options.output, 'split_{}_{}.csv'.format(i, 'test'))
        train_domains = domains[domains['contig_id'].apply(lambda contig_id: contig_id in train_contigs)]
        test_domains = domains[domains['contig_id'].apply(lambda contig_id: contig_id in test_contigs)]
        train_domains.to_csv(split_train_path, index=False)
        test_domains.to_csv(split_test_path, index=False)

    meta_path = os.path.join(options.output, 'splits.csv')
    print('Saving splits meta file to: {}'.format(meta_path))
    pd.DataFrame(splits).to_csv(meta_path, index=False)

    print('Done.')