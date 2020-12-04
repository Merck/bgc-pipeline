#!/usr/bin/env python
# David Prihoda
# Create FASTA file and Domain CSV file for positive samples used in Cimermancic et al.

import argparse
import pandas as pd
from Bio import SeqIO


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Input ClusterFinder positive domains csv file name.", metavar="FILE")
    parser.add_argument("-g", "--genome", dest="genome", required=True,
                        help="BGC genome file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name (without prefix).", metavar="FILE")
    options = parser.parse_args()

    domains = pd.read_csv(options.input)
    domains['in_cluster'] = 1
    domains.to_csv(options.output + '.csv', index=False)

    print('Writing cluster genome sequences')
    with open(options.output + '.fa', 'w') as seqfile:
        cluster_records = list(SeqIO.parse(options.genome, 'genbank'))
        SeqIO.write(cluster_records, seqfile, 'fasta')

    print('Saved to', options.output)

