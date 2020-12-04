#!/usr/bin/env python
# David Prihoda
# Input/Output utilities for reading Domain CSV files

import pandas as pd

def read_domains(file, max_evalue=None, min_bitscore=None):
    """
    Read Domain CSV file into a Domain DataFrame
    :param file: Path to Domain CSV file
    :param max_evalue: Return only domains with e-value lower than given threshold (use None to skip)
    :param min_bitscore: Return only out domains with bitscore higher than given threshold (use None to skip)
    :return: Domain DataFrame filtered by given evalue and bitscore
    """
    domains = pd.read_csv(file)
    if max_evalue:
        domains = domains[domains['evalue'] < max_evalue]
    if min_bitscore:
        if 'bitscore' not in domains:
            raise AttributeError('Cannot filter on bitscore, column not present.')
        else:
            domains = domains[domains['bitscore'] > min_bitscore]
    return domains.reset_index(drop=True)


def count_y_clusters(y):
    """
    Count BGCs regions in a list of protein domain BGC states. Done by counting all consecutive 1 as one region.
    :param y: List of BGC states (0 = non-BGC, 1 = BGC), one state for each protein domain
    :return: Count of BGCs in given list of protein domain BGC states.
    """
    prev = 0
    clusters = 0
    for val in y:
        if val == 1 and prev == 0:
            clusters += 1
        prev = val
    return clusters


def domains_to_samples(domains, sample_column='contig_id', target_column=None):
    """
    Group Domain DataFrame into a list of DataFrames and list of output states based on given sample_id column
    :param domains: Domain DataFrame with multiple sequences (samples), each marked by a unique ID in the sample_column
    :param sample_column: Sample ID column name
    :param target_column: Column name to consider as output value. If provided, will return tuple of (list of DataFrames, list of states Series).
    :return: list of DataFrames or tuple (list of DataFrames, list of states)
    """
    samples = domains.groupby(sample_column)
    sample_ids = list(samples.groups.keys())
    samples = [samples.get_group(i) for i in sample_ids]
    if target_column:
        targets = [sample[target_column] for sample in samples]
        return samples, targets
    return samples


def plot_probs(probs, window=1, ylim=(-0.1,1.1), ax=None, title='', lw=1, **kwargs):
    if not isinstance(probs, pd.Series) and not isinstance(probs, pd.DataFrame):
        probs = pd.Series(probs)
    if window !=1:
        title += ' (Rolling mean window {})'.format(window)
        probs = probs.rolling(window).mean()
    return probs.plot(figsize=None if ax else (20,2), ylim=ylim, lw=lw, ax=ax, title=title, **kwargs)