#!/usr/bin/env python
# David Prihoda
# Create filtered protein FASTA and Domain CSV files for BGC samples from MIBiG database
import argparse
import pandas as pd
import os
from Bio import SeqIO


def filter_domains(domains, mindomains, minproteins, maxevalue):
    old_size = len(domains)
    domains = domains[domains['evalue'] <= maxevalue]
    print(
        'Removed {} domains with e-value > {}, {} out of {} remained'.format(old_size - len(domains), maxevalue,
                                                                             len(domains), old_size))

    clusters = domains[['protein_id', 'contig_id']].groupby('contig_id').describe()
    num_proteins = clusters['protein_id']['unique']
    num_domains = clusters['protein_id']['count']

    old_size = len(clusters)
    small_cluster_ids = set(num_domains[num_domains < mindomains].index)
    small_cluster_ids.update(num_proteins[num_proteins < minproteins].index)
    domains = domains[domains['contig_id'].apply(lambda contig_id: contig_id not in small_cluster_ids)]
    new_size = len(domains['contig_id'].unique())
    if new_size != old_size - len(small_cluster_ids):
        raise Exception('Unexpected size {}'.format(old_size))
    print('Removed {} clusters with domains < {} or proteins < {}, {} out of {} remained'.format(old_size - new_size,
                                                                                                 mindomains,
                                                                                                 minproteins,
                                                                                                 new_size, old_size))
    return domains


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Input cluster domain csv file name.", metavar="FILE")
    parser.add_argument("-g", "--genome", dest="genome", required=True,
                        help="MIBiG genome folder path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name (without prefix).", metavar="FILE")
    parser.add_argument("-md", "--mindomains", dest="mindomains", required=True, type=int,
                        help="Minimum number of domains in cluster.", metavar="INT")
    parser.add_argument("-mp", "--minproteins", dest="minproteins", required=True, type=int,
                        help="Minimum number of proteins in cluster.", metavar="INT")
    parser.add_argument("-e", "--maxevalue", dest="maxevalue", required=True, type=float,
                        help="Maximum domain independent evalue.", metavar="FLOAT")
    parser.add_argument("-f", "--first", dest="first", action='store_true',
                        help="Only include first BGC version from each genbank file (ID suffix .1)")

    options = parser.parse_args()

    domains = pd.read_csv(options.input)

    print('Filtering clusters')
    domains = filter_domains(domains, mindomains=options.mindomains, minproteins=options.minproteins, maxevalue=options.maxevalue)

    if options.first:
        domains = domains[domains['contig_id'].apply(lambda contig_id: contig_id.endswith('.1'))]

    domains['in_cluster'] = 1
    domains.to_csv(options.output + '.csv', index=False)

    cluster_ids = domains['contig_id'].unique()

    print('Writing cluster genome sequences')
    with open(options.output + '.fa', 'w') as seqfile:
        for cluster_id_full in cluster_ids:
            cluster_id, index = cluster_id_full.split('.')
            genome_path = os.path.join(options.genome, cluster_id+'.gbk')
            cluster_records = [rec for rec in SeqIO.parse(genome_path, 'genbank') if rec.id == cluster_id_full]
            if not cluster_records:
                raise AttributeError('Cluster '+cluster_id_full+' not present in '+genome_path)
            SeqIO.write(cluster_records[0], seqfile, 'fasta')

    print('Saved to', options.output)

