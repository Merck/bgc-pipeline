#!/usr/bin/env python
# David Prihoda
# Extract the sequence and features of given BGC candidates from their annotated source bacteria GenBank files

import argparse
import pandas as pd
from Bio import SeqIO
import os
import bz2
import time

def get_candidate_genbank(cands, bacteria_path):
    """
    Create a SeqIO record for each candidate, extracted from the bzipped GenBank files in a given directory
    :param cands: DataFrame of Candidates to extract records for
    :param bacteria_path: Path to folder with bzipped GenBank bacteria files, named as contig_id.gbk.bz2
    :return: Generator of SeqIO records for each candidate
    """
    for i, candidate in cands.iterrows():
        path = os.path.join(bacteria_path, candidate['contig_id']+'.gbk.bz2')
        with bz2.open(path, 'rt') as f:
            genome = SeqIO.read(f, 'genbank')
        record = genome[candidate['nucl_start']:candidate['nucl_end']]
        # TODO assign proper locus (shorter than ~20 characters)
        record.name = 'unknown_{}'.format(i)
        record.id = candidate['candidate_id']
        record.annotations = {
            'source': '{}|{}'.format(candidate['contig_id'], candidate['species']),
            'date': time.strftime("%d-%b-%Y").upper(),
            'keywords': [candidate['candidate_hash']]
        }
        record.description = ''
        print('Yielding candidate {} with {} features'.format(i, len(record.features)))
        print(record)
        yield record


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Candidate csv file path.", metavar="FILE")
    parser.add_argument("-b", "--bacteria", dest="bacteria", required=True,
                      help="Bacteria genbank folder path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    cands = pd.read_csv(options.input)
    records = get_candidate_genbank(cands, bacteria_path=options.bacteria)

    SeqIO.write(records, options.output, 'genbank')
    print('Saved candidates to: {}'.format(options.output))

