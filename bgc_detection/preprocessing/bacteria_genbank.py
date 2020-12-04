#!/usr/bin/env python
# David Prihoda
# Create fully annotated GenBank file from bacteria Domain CSV and nucleotide fasta files.

import argparse
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
import multiprocessing
import os
import glob
import bz2

def get_bacteria_ids(bacteria_domains_path):
    path_pattern = os.path.join(bacteria_domains_path, '*.domains.csv')
    pattern_suffix = path_pattern.split('*')[-1]
    return [os.path.basename(file.replace(pattern_suffix, '')) for file in glob.glob(path_pattern)]

def run_task(args):
    i, bacteria_id, domains_path, genome_path, output_path, pfam_descriptions, evalue = args

    domains = pd.read_csv(domains_path)
    domains = domains[domains['evalue'] <= evalue]

    record = SeqIO.read(genome_path, 'fasta', alphabet=DNAAlphabet())
    record.description = record.id
    record.id = bacteria_id
    record.name = bacteria_id

    skipped = 0
    for protein_id, protein_domains in domains.groupby('protein_id', sort=False):
        protein = protein_domains.iloc[0][['gene_start', 'gene_end', 'gene_strand']]
        start = int(protein['gene_start'])
        end = int(protein['gene_end'])
        loc = FeatureLocation(start, end, strand=int(protein['gene_strand']))
        feature = SeqFeature(location=loc, qualifiers={'protein_id':[protein_id]}, type="CDS")
        record.features.append(feature)
        for d, domain in protein_domains.iterrows():
            if loc.strand == 1:
                domain_loc = FeatureLocation(start + domain['domain_start'] * 3, start + domain['domain_end'] * 3,
                                         strand=loc.strand)
            elif loc.strand == -1:
                domain_loc = FeatureLocation(end - domain['domain_end'] * 3, end - domain['domain_start'] * 3,
                                         strand=loc.strand)
            else:
                raise ValueError("Invalid strand {}".format(loc.strand))
            if domain_loc.start < loc.start or domain_loc.start > loc.end or domain_loc.end < loc.start or domain_loc.end > loc.end:
                skipped += 1
                #print("WARNING: Skipping invalid domain {} with location: {} outside protein location {}.".format(domain['pfam_id'], domain_loc, loc))
                continue
            domain_feature = SeqFeature(location=domain_loc, type="PFAM_domain",
                                        qualifiers={
                                            'pfam_id': domain['pfam_id'],
                                            'evalue': domain['evalue'],
                                            'description': [pfam_descriptions.loc[domain['pfam_id']]]})
            record.features.append(domain_feature)
    if skipped:
        print("WARNING: Skipped {} domains with invalid location outside protein location".format(skipped))

    with bz2.open(output_path, 'wt') as f:
        SeqIO.write(record, f, 'genbank')

    print('Saved {} features (task_{}) to: {}'.format(len(record.features), i, output_path))

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input",
                      help="Input bacteria domain csv folder path.", metavar="FILE")
    parser.add_argument("-g", "--genome", dest="genome",
                      help="Input bacteria genome fasta folder path.", metavar="FILE")
    parser.add_argument("-d", "--description", dest="description",
                      help="Pfam descriptions csv file.", metavar="FILE")
    parser.add_argument("-e", "--maxevalue", dest="maxevalue", required=True, type=float,
                        help="Maximum domain independent e-value.", metavar="FLOAT")
    parser.add_argument("-o", "--output", dest="output",
                      help="Output genbank folder name.", metavar="FILE")
    options = parser.parse_args()

    bacteria_ids = get_bacteria_ids(options.input)
    pfam_descriptions = pd.read_csv(options.description).set_index('pfam_id')['description']

    print('Creating output directory: {}'.format(options.output))
    os.makedirs(options.output)

    tasks = [(
                  i,
                  bacteria_id,
                  os.path.join(options.input, bacteria_id + '.domains.csv'),
                  os.path.join(options.genome, bacteria_id + '.genome.fa'),
                  os.path.join(options.output, bacteria_id + '.gbk.bz2'),
                  pfam_descriptions,
                  options.maxevalue
              ) for i, bacteria_id in enumerate(bacteria_ids)]

    pool = multiprocessing.Pool()
    pool.map(run_task, tasks)
    print('Done.')

