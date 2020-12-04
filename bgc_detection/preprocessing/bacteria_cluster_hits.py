#!/usr/bin/env python
# David Prihoda
# Convert QcHitsRanno.tsv into ClusterHits.csv with simplified columns

import pandas as pd
import argparse

def get_bacteria_id(genome_id):
    return genome_id.split('|')[2].split('_')[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input",
                      help="Input tsv file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output",
                      help="Output csv file path.", metavar="FILE")
    options = parser.parse_args()

    clusters = pd.read_csv(options.input, sep='\s+')

    # Get short bacteria ID
    clusters['bacteria_id'] = clusters['Genome_ID'].map(lambda genome_id: get_bacteria_id(genome_id))
    clusters.drop('Genome_ID', axis=1, inplace=True)

    # Get strand based on order of start-end
    clusters['strand'] = (clusters['Genome_Start'] < clusters['Genome_End']).map(lambda b: 1 if b else -1)

    # Flip start and end in -1 strand
    clusters['cluster_start'] = clusters.apply(lambda row: row['Genome_Start'] if row['strand'] == 1 else row['Genome_End'], axis=1)
    clusters['cluster_end'] = clusters.apply(lambda row: row['Genome_End'] if row['strand'] == 1 else row['Genome_Start'], axis=1)
    clusters['BGC_ID'] = clusters['BGC_ID'].apply(lambda bgc_id: bgc_id.split('_')[0])
    clusters['contig_id'] = clusters['BGC_ID'].apply(lambda bgc_id: '.'.join(bgc_id.split('_')[:2]))

    clusters = clusters[['BGC_ID','bacteria_id','strand','cluster_start','cluster_end']]

    clusters.to_csv(options.output, index=False)