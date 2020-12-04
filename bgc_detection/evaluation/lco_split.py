#!/usr/bin/env python
# David Prihoda
# Create Leave-Class-Out splits from set of negative and positive samples
# Produces a folder with (train, test) files for each class.

import merge_split_samples
import argparse
import numpy as np
import pandas as pd
import os
from sklearn.model_selection import KFold, ShuffleSplit

NEG_CLASS_NAME = '_neg_'

class LeaveClassOutSplitter:
    """
    Splitter that splits Series of samples in Leave-Class-Out fashion.
    """

    def __init__(self, unique_classes, neg_class: str, neg_test_size, pos_test_count=None, random_state=0):
        """

        :param unique_classes: Unique positive classes that will be present when splitting
        :param neg_class: Class label to handle as negative class
        :param neg_test_size: Fraction of negative samples to use for testing. The rest is used for training.
        :param pos_test_count: Upsample positive samples to given number. If pos_test_count is lower than the number
        of samples in given class, sampling without replacement is used. If it is higher, we first select all of the samples of given class
        and sample the remainder with replacement.
        :param random_state: Random state to use for sampling.
        """
        self.unique_classes = unique_classes
        self.neg_class = neg_class
        if neg_test_size < 0 or neg_test_size > 1:
            raise AttributeError('Negative test size has to be a fraction between 0.0 and 1.0')
        self.neg_splitter = ShuffleSplit(test_size=neg_test_size, random_state=random_state)
        self.pos_test_count = pos_test_count

    def split(self, samples, classes):
        """
        Split Series of samples in Leave-Class-Out fashion.
        Returns a generator that produces a (train, test) tuple for each positive class.
        The train split contains a random subset of negative samples (marked as self.neg_class) and all positive samples except given class.
        The test split contains the remaining samples from the negative set and all samples from given class.
        Positive samples can be upsampled by using pos_test_count. If pos_test_count is lower than the number
        of samples in given class, sampling without replacement is used. If it is higher, we first select all of the samples of given class
        and sample the remainder with replacement.

        :param samples: Pandas Series of DataFrames (samples)
        :param classes: Series or list of classes for each sample
        :return: Generator of (train, test) split indexes, will generate one tuple for each positive class
        """
        if len(samples) != len(classes):
            raise AttributeError("Samples and classes have to be the same length")
        neg_idx = np.where(classes == self.neg_class)[0]
        if not len(neg_idx):
            raise AttributeError("No negative samples. Add samples with class = {}.".format(self.neg_class))

        neg_samples = samples[neg_idx]
        neg_train_idx, neg_test_idx = next(self.neg_splitter.split(neg_samples))
        neg_train_idx = neg_idx[neg_train_idx]
        neg_test_idx = neg_idx[neg_test_idx]

        for klass in self.unique_classes:
            # Train on all other classes except negative
            pos_train_idx = np.where((classes != klass) & (classes != self.neg_class))[0]
            # Test on given class
            pos_test_idx = np.where(classes == klass)[0]
            num_class_samples = len(pos_test_idx)
            if not num_class_samples:
                print('No samples of class {} found'.format(klass))
            if self.pos_test_count:
                if num_class_samples > self.pos_test_count:
                    # we have more samples than we are sampling, choose without replacement
                    pos_test_idx = np.random.choice(pos_test_idx, self.pos_test_count, replace=False)
                else:
                    # we have less samples than we are sampling, use all indexes + a sampled remainder
                    num_remaining = self.pos_test_count - num_class_samples
                    remaining_choice = np.random.choice(pos_test_idx, num_remaining, replace=True)
                    pos_test_idx = np.concatenate([pos_test_idx, remaining_choice])
            print('Train: {} pos, {} neg. Test: {} pos, {} neg'.format(len(pos_train_idx), len(neg_train_idx), len(pos_test_idx), len(neg_test_idx)))
            # Return unions of negative and positive splits
            yield np.concatenate([pos_train_idx, neg_train_idx]), np.concatenate([pos_test_idx, neg_test_idx])


def filter_lco_samples(pos_samples, pos_classes, neg_samples):
    # Remove hybrids, Other and unknown
    selected_classes = sorted(list(set([c for c in set(pos_classes) if ';' not in c and 'Other' not in c and '?' not in c and 'Nucleoside' not in c])))
    print('Selected classes:', selected_classes)

    # Get positive samples of selected classes
    sel_pos_samples = [s for s, c in zip(pos_samples, pos_classes) if c in selected_classes]
    sel_pos_classes = [c for c in pos_classes if c in selected_classes]
    print('{} non-hybrid of {} total BGCs remained'.format(len(sel_pos_samples), len(pos_samples)))

    # Merge into one list
    lco_samples = sel_pos_samples + neg_samples
    lco_classes = np.concatenate([sel_pos_classes, np.array([NEG_CLASS_NAME] * len(neg_samples))])
    print('{} total samples'.format(len(lco_samples)))

    return selected_classes, lco_samples, lco_classes


def lco_samples(selected_classes, val_samples, val_classes, result_path, random_seeds, pos_test_count=300):
    meta = []
    for random_seed in random_seeds:
        print('Random seed', random_seed)
        lco_splitter = LeaveClassOutSplitter(
            unique_classes=selected_classes,
            neg_class=NEG_CLASS_NAME,
            neg_test_size=0.33,
            pos_test_count=pos_test_count,
            random_state=random_seed
        )
        labels = ['{}\n({} BGCs sampled {}x)'.format(class_name, sum(val_classes == class_name), pos_test_count) for class_name in selected_classes]

        np.random.seed(random_seed)
        merged_splits = merge_split_samples.merged_split(
            val_samples,
            lco_splitter,
            shuffle_train=False,
            shuffle_test=True,
            split_params={'classes': val_classes}
        )
        merged_splits = list(merged_splits)

        for split_domains, class_name, label in zip(merged_splits, selected_classes, labels):
            split_name = '{}.seed{}'.format(class_name, random_seed)
            train_domains, test_domains = split_domains
            train_csv_path = os.path.join(result_path, split_name+'.train.csv')
            train_domains.to_csv(train_csv_path, index=False)
            print('Saved LCO {} train sequence to: {}'.format(split_name, train_csv_path))
            test_csv_path = os.path.join(result_path, split_name+'.test.csv')
            test_domains.to_csv(test_csv_path, index=False)
            print('Saved LCO {} test sequence to: {}'.format(split_name, test_csv_path))
            meta.append({
                'label': label,
                'name': split_name,
                'group': class_name
            })

    meta_csv_path = os.path.join(result_path, 'splits.csv')
    meta = pd.DataFrame(meta)
    meta.to_csv(meta_csv_path, index=False)
    print('Saved splits meta file to: {}'.format(meta_csv_path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--maxevalue", dest="maxevalue", required=True, type=float,
                        help="Maximum domain independent e-value.", metavar="FLOAT")
    parser.add_argument("-p", "--positive", dest="positive", required=True,
                        help="Path to positive samples file.", metavar="FILE")
    parser.add_argument("-c", "--classes", dest="classes", required=True,
                        help="Path to csv file containing classes of positive samples.", metavar="FILE")
    parser.add_argument("--classes-column", dest="classes_column", default="classes",
                        help="Class column in classes file.", metavar="STRING")
    parser.add_argument("--pos-test-count", dest="pos_test_count", required=False, default=300, type=int,
                        help="Number of positive test samples (use sampling with replacement).", metavar="INT")
    parser.add_argument("--random-seed", dest="random_seed", required=False, default=[], type=int, action='append',
                        help="Random seed used to shuffle the samples.", metavar="INT")
    parser.add_argument("-n", "--negative", dest="negative", required=True,
                        help="Path to negative samples file.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output samples folder path.", metavar="FILE")
    options = parser.parse_args()

    pos_domains = pd.read_csv(options.positive)
    pos_domains = pos_domains[pos_domains['evalue'] < options.maxevalue]
    pos_samples = [s for i, s in pos_domains.groupby('contig_id')]
    pos_ids = np.array([sample['contig_id'].iloc[0] for sample in pos_samples])

    neg_domains = pd.read_csv(options.negative)
    neg_domains = neg_domains[neg_domains['evalue'] < options.maxevalue]
    neg_samples = [s for i, s in neg_domains.groupby('contig_id')]

    properties = pd.read_csv(options.classes).set_index('contig_id')
    pos_classes = properties[options.classes_column][pos_ids]

    selected_classes, val_samples, val_classes = filter_lco_samples(
        pos_samples=pos_samples,
        pos_classes=pos_classes,
        neg_samples=neg_samples
    )

    print('Output will be saved to {}/'.format(os.path.abspath(options.output)))
    os.makedirs(options.output)

    if not options.random_seed:
        options.random_seed = [0]

    lco_samples(
        selected_classes=selected_classes,
        val_samples=val_samples,
        val_classes=val_classes,
        result_path=options.output,
        pos_test_count=options.pos_test_count,
        random_seeds=options.random_seed
    )
