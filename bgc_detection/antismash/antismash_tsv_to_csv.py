#!/usr/bin/env python
# David Prihoda
# Convert an antiSMASH TSV file taken from 'txt/*_BGC.txt' file into a Candidate CSV file

import argparse
import pandas as pd
import numpy as np

def antismash_tsv_candidates(export):
    """
    Convert an antiSMASH TSV file taken from 'txt/*_BGC.txt' into a Candidate DataFrame
    :param export: DataFrame of antiSMASH TSV file taken from 'txt/*_BGC.txt'
    :return: Candidate DataFrame
    """
    candidates = pd.DataFrame()
    candidates['contig_id'] = export['Contig']
    candidates['nucl_start'] = export['BGC_range'].apply(lambda r: r.split(';')[0]).astype(np.int)
    candidates['nucl_end'] = export['BGC_range'].apply(lambda r: r.split(';')[1]).astype(np.int)
    candidates['nucl_length'] = candidates['nucl_end'] - candidates['nucl_start'] + 1
    candidates['candidate_id'] = candidates.apply(lambda row: '{}:{}-{}'.format(row['contig_id'], row['nucl_start'], row['nucl_end']), axis=1)
    candidates['classes'] = export['BGC type']
    candidates['detection_rules'] = export['detection rules used']
    candidates['genes'] = export['genes']
    candidates['subclusters'] = export['subclusters']
    candidates['NRPSs/PKSs'] = export['NRPSs/PKSs']
    candidates['signature_genes'] = export['signature_genes']
    candidates['RiPPs'] = export['RiPPs']
    candidates['predicted_structure'] = export['predicted structure']
    candidates['monomers'] = export['monomers']
    return candidates

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Path to AntiSMASH txt result (merged into tsv) file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    export = pd.read_csv(options.input, sep='\t')
    candidates = antismash_tsv_candidates(export)
    candidates.to_csv(options.output, index=False)

    print('Saved {} candidates to {}'.format(len(candidates), options.output))

