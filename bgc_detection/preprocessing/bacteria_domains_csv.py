#!/usr/bin/env python
# David Prihoda
# Create a Domain CSV file from bacteria hmmscan domtbl file and the protein fasta file
# Protein fasta file is needed to obtain gene coordinates which are not present in hmmscan file

import argparse
from Bio import SeqIO
from domtbl2csv import domtbl_to_df, normalize_gene_coord
import multiprocessing
import os
import glob

def get_location_from_description(descr):
    s = descr.split('#')
    # count from the end to avoid errors when '#' is present in id
    return int(s[-4]), int(s[-3]), int(s[-2])

def get_gene_locations(protein_path):
    # Note: prodigal fasta file needs to be used (location is included between hashtags in description)
    proteins = SeqIO.parse(protein_path, 'fasta')
    gene_locations = {record.id : get_location_from_description(record.description) for record in proteins}

    return gene_locations


def get_bacteria_domains_df(domains, gene_locations):
    domains['protein_id'] = domains['sequence_id'].apply(lambda s: s.split('|')[1]+'_'+s.split('_')[-1])
    domains['gene_start'] = domains['sequence_id'].apply(lambda g: normalize_gene_coord(gene_locations[g][0]))
    domains['gene_end'] = domains['sequence_id'].apply(lambda g: normalize_gene_coord(gene_locations[g][1]))
    domains['gene_strand'] = domains['sequence_id'].apply(lambda g: gene_locations[g][2])
    return domains

def get_bacteria_ids(bacteria_domtbl_path):
    path_pattern = os.path.join(bacteria_domtbl_path, '*.domtbl.tbl')
    pattern_suffix = path_pattern.split('*')[-1]
    return [os.path.basename(file.replace(pattern_suffix, '')) for file in glob.glob(path_pattern)]

def run_task(args):
    i, domtbl_path, protein_path, output_path = args

    domains = domtbl_to_df(domtbl_path)
    gene_locations = get_gene_locations(protein_path)
    df = get_bacteria_domains_df(domains, gene_locations)

    df = df[['protein_id', 'gene_start', 'gene_end', 'gene_strand',
             'pfam_id', 'domain_start', 'domain_end', 'evalue', 'bitscore']]
    df.to_csv(output_path, index=False)

    print('Saved task_{} to: {}'.format(i, output_path))

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input",
                      help="Input bacteria hmmscan domtbl folder.", metavar="FILE")
    parser.add_argument("-p", "--protein", dest="protein",
                      help="Input bacteria protein fasta folder.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output",
                      help="Output csv folder name.", metavar="FILE")
    options = parser.parse_args()

    bacteria_ids = get_bacteria_ids(options.input)

    tasks = [(
                  i,
                  os.path.join(options.input, bacteria_id+'.domtbl.tbl'),
                  os.path.join(options.protein, bacteria_id + '.proteins.fa'),
                  os.path.join(options.output, bacteria_id+'.domains.csv')
              ) for i, bacteria_id in enumerate(bacteria_ids)]

    pool = multiprocessing.Pool()
    pool.map(run_task, tasks)
    print('Done.')

