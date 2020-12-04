#!/usr/bin/env python
# David Prihoda
# Average predictions in a Domain CSV file by protein

import argparse
import pandas as pd

PROTEIN_GROUP_COLS = ['gene_start', 'gene_end', 'gene_strand', 'protein_id']
PROTEIN_EXTRA_COLS = ['contig_id', 'in_cluster']

def agg_concat(s):
    return ';'.join(s)

def average_protein_prediction(domains, y=None, concat_domains=True):
    """
    Average predictions into a 'prediction' column by protein using the 'protein_id' and other PROTEIN_GROUP_COLS.
    :param domains: DataFrame from the Domain CSV file
    :param y: Series of predictions to be averaged and written in the 'prediction' column
    :param concat_domains: Whether to include a ';'-concatenated list of pfam_ids for each protein.
    :return: DataFrame of proteins with averaged 'prediction' column
    """
    extra_cols = [col for col in PROTEIN_EXTRA_COLS if col in domains.columns]
    cols = extra_cols + PROTEIN_GROUP_COLS
    if concat_domains:
        cols.append('pfam_id')
    copy = domains[cols].copy()
    copy['prediction'] = y
    per_gene = copy.groupby(extra_cols + PROTEIN_GROUP_COLS, sort=False)
    if concat_domains:
        return per_gene.agg({'pfam_id': agg_concat, 'prediction': 'mean'})\
            .rename(columns={'pfam_id': 'pfam_ids'})\
            .reset_index()
    else:
        return per_gene.mean().reset_index()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Path to domain prediction CSV.", metavar="FILE")
    parser.add_argument("-c", "--column", dest="column", default='prediction',
                        help="Prediction column.", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file path.", metavar="FILE")
    options = parser.parse_args()

    domains = pd.read_csv(options.input)

    proteins = average_protein_prediction(domains, domains[options.column])

    proteins.to_csv(options.output, index=False)
    print('Saved protein predictions to: {}'.format(options.output))