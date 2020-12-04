#!/usr/bin/env python
# David Prihoda
# Prepare labelled Domain CSV file for the Cimermancic et al. annotated contigs
# Uses the hmmscan domtbl file and the CSV containing positive (BGC) genes

import argparse
from domtbl2csv import domtbl_to_df
import pandas as pd


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input",
                      help="Input cluster hmmscan domtbl file path.", metavar="FILE")
    parser.add_argument("-l", "--labels", dest="labels",
                      help="BGC gene CSV file path with 'contig_id' and 'locus_tag' columns. "
                           "Used to mark positive protein domains.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output",
                      help="Output csv file path.", metavar="FILE")
    options = parser.parse_args()

    df = domtbl_to_df(options.input, format='proteins2fasta')

    true_cluster_genes = pd.read_csv(options.labels, sep=';')
    true_cluster_genes['in_cluster'] = 1

    df = df.merge(true_cluster_genes[['contig_id','locus_tag','in_cluster']], on=['contig_id','locus_tag'], how='left').fillna(0)

    df[['contig_id','locus_tag','protein_id','gene_start','gene_end','gene_strand','pfam_id','domain_start','domain_end','evalue','bitscore','in_cluster']]\
        .to_csv(options.output, index=False)

    print('Saved to', options.output)

