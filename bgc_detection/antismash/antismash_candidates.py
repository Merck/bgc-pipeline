#!/usr/bin/env python
# David Prihoda
# Convert the final antiSMASH GenBank file to a BGC Candidate CSV file.

import argparse
import pandas as pd
from Bio import SeqIO

def antismash_candidates(records):
    """
    Create a DataFrame from SeqIO records from the final antiSMASH GenBank file.
    :param records: SeqIO records from the final antiSMASH GenBank file.
    :return: DataFrame with BGC candidates
    """
    antismash = []
    for record in records:
        clusters = [f for f in record.features if f.type == 'cluster']
        antismash += [{
            'nucl_start': int(cluster.location.start),
            'nucl_end': int(cluster.location.end),
            'classes': cluster.qualifiers.get('product', [None])[0],
            'contig_id': record.id
        } for cluster in clusters]
    return pd.DataFrame(antismash)

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Path to AntiSMASH final GenBank file.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    records = SeqIO.parse(options.input, 'genbank')

    candidates = antismash_candidates(records)
    candidates.to_csv(options.output, index=False)

    print('Saved {} candidates to {}'.format(len(candidates), options.output))

