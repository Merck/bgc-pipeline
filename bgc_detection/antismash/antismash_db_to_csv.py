#!/usr/bin/env python
# David Prihoda
# Convert TSV result downloaded using antismash_db_export.py script into a Candidate CSV file

import argparse
import pandas as pd

EXPORT_COLUMNS = ['Genus', 'Species', 'NCBI accession', 'Cluster number', 'BGC type', 'From', 'To', 'Most similar known cluster','Similarity in %','MIBiG BGC-ID','Results URL']

def read_antismash_export(path):
    df = pd.read_csv(path, sep='\t', comment='#', header=None)
    df.columns = EXPORT_COLUMNS
    return df

def antismash_export_candidates(export):
    """
    Convert TSV result downloaded using antismash_db_export.py script into a Candidate CSV file
    :param export: DataFrame with TSV result from antiSMASH database
    :return: Candidate DataFrame
    """
    candidates = pd.DataFrame()
    candidates['contig_id'] = export['NCBI accession']
    candidates['candidate_id'] = export.apply(lambda row: '{}:{}'.format(row['NCBI accession'], row['Cluster number']), axis=1)
    candidates['nucl_start'] = export['From']
    candidates['nucl_end'] = export['To']
    candidates['classes'] = export['BGC type']
    candidates['similar_bgc_id'] = export['MIBiG BGC-ID'].apply(lambda x: x.replace('_c', '.') if x and x != 'None' else None)
    candidates['similar_percent'] = export['Similarity in %'].apply(lambda x: float(x) if x and x != 'None' else None)
    return candidates

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Path to AntiSMASH DB export TSV file.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    export = read_antismash_export(options.input)
    candidates = antismash_export_candidates(export)
    candidates.to_csv(options.output, index=False)

    print('Saved {} candidates to {}'.format(len(candidates), options.output))

