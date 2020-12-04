#!/usr/bin/env python
# David Prihoda
# Functions for shuffling and merging sequence samples (DataFrames) into a long sequence
# Used for fake genome generation from positive and negative BGC samples

import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

def merge_samples(samples_series, idx, shuffle=True):
    """
    Merge a Series of DataFrames into a single DataFrame, filtered by a Series index. DataFrames can be shuffled before merging.
    Used to generate an artificial genome (long DataFrame) from a series of samples (short DataFrames)

    :param samples_series: Series of DataFrames
    :param idx: Array of indexes used to select DataFrames for merging
    :param shuffle: Whether to shuffle the DataFrames
    :return: Subset of given series selected by idx, shuffled if specified and merged to a single DataFrame
    """
    if shuffle:
        np.random.shuffle(idx)
    return pd.concat(list(samples_series[idx]))


def merged_split(samples_list, splitter, shuffle_train=True, shuffle_test=True, split_params=None):
    """
    Create generator of random train and test splits, where each train and test split
    is a single DataFrame created from shuffled and merged samples using merge_samples function.
    Will generate given number of (train, test) splits based on splitter argument.

    :param samples_list: list of DataFrames to repeatedly split and merge
    :param splitter: Number of KFold splits or Splitter with split(samples) function that will be used
    :param shuffle_train: Whether to shuffle the samples in the train split before merging
    :param shuffle_test: Whether to shuffle the samples in the test split before merging
    :param split_params: Additional arguments to pass to the splitter split function
    :return: Generator of (train, test) splits for given list of samples, where each train and test split
    is a single DataFrame created from shuffled and merged samples using merge_samples function.
    """
    if split_params is None:
        split_params = {}
    if isinstance(splitter, int):
        splitter = KFold(n_splits=splitter, shuffle=True)
    indexable_X = pd.Series(samples_list)
    for trainidx, testidx in splitter.split(indexable_X, **split_params):
        train_X = merge_samples(indexable_X, trainidx, shuffle=shuffle_train)
        test_X = merge_samples(indexable_X, testidx, shuffle=shuffle_test)
        yield train_X, test_X
