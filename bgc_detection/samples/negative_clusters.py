#!/usr/bin/env python
# David Prihoda
# Generate artificial negative training set from a set of reference bacteria and a set of positive samples

import argparse
import pandas as pd
import numpy as np
from positive_mibig import filter_domains
from os.path import basename
import glob
import os
from multiprocessing import Pool

def get_gene_summary(protein_id, sequence_group):
    first = sequence_group.iloc[0]
    return {
        # todo: NOTE: if gene domains are counted in the future, bacteria domains have to be filtered by evalue first
        'protein_id': protein_id,
        'num_domains': len(sequence_group),
        'nucl_length': abs(first['gene_end'] - first['gene_start']) + 1,
        'gene_start': first['gene_start'],
        'gene_end': first['gene_end'],
        'gene_strand': first.get('gene_strand')
    }


def get_gene_summaries(domains):
    return pd.DataFrame([get_gene_summary(i, group) for i, group in domains.groupby('protein_id')])


def get_similar_percent(value, sorted_df, key, similar_fraction):
    size = len(sorted_df)
    interval_half_with = max(1, round(similar_fraction * size / 2))
    num_smaller = sum(sorted_df[key] <= value)
    left_index = max(0, num_smaller - interval_half_with)
    right_index = min(size, num_smaller + interval_half_with)
    return sorted_df[left_index:right_index]


def check_overlap(first_start, first_end, second_start, second_end):
    astart = min([first_start, first_end])
    aend = max([first_start, first_end])
    bstart = min([second_start, second_end])
    bend = max([second_start, second_end])
    no_overlap = aend <= bstart or astart >= bend
    return not no_overlap


def get_bacteria_ids(bacteria_domains_path):
    path_pattern = os.path.join(bacteria_domains_path, '*.domains.csv')
    pattern_suffix = path_pattern.split('*')[-1] # suffix such as '.domains.csv' will be removed to get bacteria id
    return [basename(file.replace(pattern_suffix, '')) for file in glob.glob(path_pattern)]


def create_negative_clusters(i, bgc_domains, bacteria_id, bacteria_domains, all_hits, seed, num_samples, similar_fraction):
    """
    Generate given number of negative samples for given reference bacteria. Each negative sample will be based on a random positive sample from given set.
    Each gene in the positive sample is replaced with a random gene from the reference bacteria,
    while considering only given fraction of genes that were most similar in number of Pfam domains.
    Genes within given coordinates (hits) are removed.
    :param i: Index used for logging purposes
    :param bgc_domains: Merged DataFrame with all positive samples, each sample is distinguished using a unique 'contig_id' value
    :param bacteria_id: Id of reference bacteria to get genes from
    :param bacteria_domains: DataFrame of domains in reference bacteria
    :param all_hits: Coordinates to remove from reference bacteria before sampling random genes.
    :param seed: Random seed.
    :param num_samples: Number of samples to generate
    :param similar_fraction: Sample from given fraction of most similar genes by number of Pfam domains. Genes are sorted by number of domains and then only similar_fraction / 2 above and below are considered.
    :return: DataFrame of negative samples, each sample has a unique 'contig_id' value prefixed with 'NEG_FAKE_CLUSTER|'
    """
    sample_ids = set()

    bgc_ids = bgc_domains['contig_id'].unique()
    domains_by_contig_id = bgc_domains.groupby('contig_id')
    num_bacteria = len(all_bacteria_ids)

    # Seed using array of numbers based on global seed and each char of bacteria ID
    bac_seed = [seed] + [ord(c) for c in bacteria_id]
    rand = np.random.RandomState(bac_seed)
    random_bgc_ids = rand.choice(bgc_ids, size=num_samples, replace=False)
    print('Generating {} samples for bacteria {} ({}/{}), using BGCs: {}'.format(num_samples, bacteria_id, i, num_bacteria, random_bgc_ids))

    genepool = get_gene_summaries(bacteria_domains)
    genepool_by_num_domains = genepool.sort_values(by='num_domains')

    if all_hits is not None:
        bacteria_hits = all_hits[all_hits['bacteria_id'].apply(lambda i: i == bacteria_id or i == bacteria_id.split('.')[0])]
        num_hits = len(bacteria_hits)
        orig_len = len(bacteria_domains)
        for hit in bacteria_hits.to_records():
            bacteria_domains = bacteria_domains[bacteria_domains.apply(
                lambda b: not check_overlap(b['gene_start'], b['gene_end'], hit['cluster_start'],
                                            hit['cluster_end']), axis=1)]
        removed = orig_len - len(bacteria_domains)
        print(' Removing {} known inter-cluster regions from bacteria {}: {}/{} removed = {:.1f}%'.format(num_hits, bacteria_id, removed, orig_len, removed / orig_len * 100))

    result_samples = []
    for full_cluster_id in random_bgc_ids:
        cluster_domains = domains_by_contig_id.get_group(full_cluster_id)
        sample_id = 'NEG_FAKE_CLUSTER|' + bacteria_id + '|' + full_cluster_id
        if sample_id in sample_ids:
            raise Exception('Duplicate sample ID found: {}'.format(sample_id))
        sample_ids.add(sample_id)

        negative_genes = []
        curr_length = 0
        for protein_id, protein_group in cluster_domains.groupby('protein_id'):
            gene = get_gene_summary(protein_id, protein_group)

            similar_genes = get_similar_percent(gene['num_domains'], genepool_by_num_domains, 'num_domains',
                                                      similar_fraction)
            replacement = similar_genes.sample(1, random_state=rand).iloc[0].copy()
            replacement_length = replacement['gene_end'] - replacement['gene_start']
            #print('  Found {} most similar genes by num_domains, replacing {} domains with {} domains'.format(
            #    len(num_domains_similar), gene['num_domains'], replacement['num_domains']))

            replacement['gene_start'] = curr_length
            replacement['gene_end'] = curr_length + replacement_length
            curr_length = replacement['gene_end']

            negative_genes.append(replacement)

        negative_genes = pd.DataFrame(negative_genes)
        negative_domains = negative_genes.merge(
            bacteria_domains[['protein_id', 'pfam_id', 'domain_start', 'domain_end', 'evalue', 'bitscore']],
            on='protein_id', how='left')
        negative_domains['contig_id'] = sample_id

        negative_domains['in_cluster'] = 0
        negative_domains = negative_domains[['contig_id', 'protein_id', 'gene_start', 'gene_end', 'gene_strand',
                                             'pfam_id', 'domain_start', 'domain_end', 'evalue', 'bitscore','in_cluster']]
        result_samples.append(negative_domains)

    return result_samples


def run_task(args):
    i, bacteria_id, bgc_domains, all_hits, options = args
    bacteria_domains = pd.read_csv(os.path.join(options.bacteria, bacteria_id + '.domains.csv'))
    bacteria_domains = bacteria_domains[bacteria_domains['evalue'] <= options.max_evalue]
    return create_negative_clusters(
        i=i,
        bgc_domains=bgc_domains,
        bacteria_id=bacteria_id,
        bacteria_domains=bacteria_domains,
        all_hits=all_hits,
        seed=options.random,
        num_samples=options.num,
        similar_fraction=options.similar_fraction
    )


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Input BGC domain CSV file name.", metavar="FILE")
    parser.add_argument("-b", "--bacteria", dest="bacteria", required=True,
                        help="Bacteria CSV domains folder path.", metavar="FILE")
    parser.add_argument("-bh", "--hits", dest="hits", required=False,
                        help="Bacteria BGC hits CSV file path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output CSV file path.", metavar="FILE")
    parser.add_argument("-md", "--min-domains", dest="min_domains", required=True, type=int,
                        help="Minimum number of domains in cluster.", metavar="INT")
    parser.add_argument("-mp", "--min-proteins", dest="min_proteins", required=True, type=int,
                        help="Minimum number of proteins in cluster.", metavar="INT")
    parser.add_argument("-e", "--max-evalue", dest="max_evalue", required=True, type=float,
                        help="Maximum domain independent evalue.", metavar="FLOAT")
    parser.add_argument("-s", "--similar-fraction", dest="similar_fraction", required=True, type=float,
                        help="Percentage of most similar genes to choose from.", metavar="FLOAT")
    parser.add_argument("-n", "--num", dest="num", required=True, type=int,
                        help="Number of samples per bacteria.", metavar="INT")
    parser.add_argument("-r", "--random", dest="random", required=False, type=int, default=123,
                        help="Random seed used for picking genes and bacteria.", metavar="INT")

    options = parser.parse_args()

    bgc_domains = pd.read_csv(options.input)
    bgc_domains = filter_domains(
        bgc_domains,
        mindomains=options.min_domains,
        minproteins=options.min_proteins,
        maxevalue=options.max_evalue
    )

    all_bacteria_ids = get_bacteria_ids(options.bacteria)
    print('Using corpus of {} bacteria'.format(len(all_bacteria_ids)))

    all_hits = pd.read_csv(options.hits) if options.hits else None

    tasks = [(i, bacteria_id, bgc_domains, all_hits, options) for i, bacteria_id in enumerate(all_bacteria_ids)]

    pool = Pool()
    sample_lists = pool.map(run_task, tasks)
    samples: pd.DataFrame = pd.concat([s for samples in sample_lists for s in samples])

    print('Saving {} negative samples with {} domains to {}'.format(len(samples['contig_id'].unique()), len(samples), options.output))
    samples.to_csv(options.output, index=False)
    print('Done.')