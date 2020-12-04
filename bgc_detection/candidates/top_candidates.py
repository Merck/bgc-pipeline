#!/usr/bin/env python
# David Prihoda
# Merge different candidate properties into a single CSV file

import argparse
import pandas as pd
import os

def get_clusterfinder_prediction(candidate, predictions_dir):
    """
    Get average ClusterFinder HMM prediction of domains in the location of the given candidate
    :param candidate: single BGC candidate Series with nucl_start, nucl_end and contig_id values
    :param predictions_dir: Path to ClusterFinder HMM predictions directory, each domain CSV file should be named contig_id.csv
    :return: Average value of ClusterFinder HMM prediction of domains in the location of the given candidate
    """
    predictions = pd.read_csv(os.path.join(predictions_dir, candidate['contig_id']+'.csv'))
    predictions = predictions[predictions['gene_start'] >= candidate['nucl_start']]
    predictions = predictions[predictions['gene_end'] <= candidate['nucl_end']]
    return predictions['prediction'].mean()

def get_top_candidates(candidates, clusterfinder_dir, species, merge):
    """
    Merge different candidate properties into a single DataFrame using the 'candidate_id' column
    :param candidates: Reference DataFrame of candidates
    :param clusterfinder_dir: Path to ClusterFinder HMM predictions directory, each domain CSV file should be named contig_id.csv
    :param species: DataFrame of species
    :param merge: List of DataFrames to merge using the 'candidate_id' column
    :return: Candidates merged into a single DataFrame
    """
    merged = None
    for i, df in enumerate(merge):
        merged = merged.join(df, how='outer', rsuffix=str(i)) if merged is not None else df
    # prepend species column to candidates table
    candidates: pd.DataFrame = species.merge(candidates, on='contig_id', how='right').set_index('candidate_hash')
    # TODO: inefficient implementation that opens the prediction file again for each candidate
    # TODO: Could be implemented by grouping by bacteria first
    cf_prediction = candidates.apply(lambda candidate: get_clusterfinder_prediction(candidate, clusterfinder_dir), axis=1)
    candidates.insert(2, 'cf_prediction', cf_prediction)
    merged = merged.join(candidates, how='inner')
    merged.index.name = 'candidate_hash'
    return merged

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Target model candidate csv file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")
    parser.add_argument("-c", "--clusterfinder", dest="clusterfinder", required=True,
                      help="ClusterFinder domain-level prediction CSV folder path.", metavar="FILE")
    parser.add_argument("-s", "--species", dest="species", required=True,
                      help="Bacterial species tsv file path.", metavar="FILE")
    parser.add_argument(dest='merge', nargs='+',
                        help="Paths to additional csv files with 'contig_id' column to merge with candidate csv file "
                             "on the 'candidate_hash' column.", metavar="FILE")

    options = parser.parse_args()

    merged_candidates = get_top_candidates(
        candidates=pd.read_csv(options.input),
        species=pd.read_csv(options.species, sep='\t'),
        clusterfinder_dir=options.clusterfinder,
        merge=[pd.read_csv(path).set_index('contig_id') for path in options.merge]
    )

    merged_candidates.to_csv(options.output, index=True)
    print('Saved {} top candidates to: {}'.format(len(merged_candidates), options.output))

