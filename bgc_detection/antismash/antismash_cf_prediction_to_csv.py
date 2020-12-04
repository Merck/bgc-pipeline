#!/usr/bin/env python
# David Prihoda
# Convert the final antiSMASH GenBank file with ClusterFinder predictions into a Domain CSV with predictions

import argparse
import pandas as pd
from Bio import SeqIO

def get_qualifier_without_prefix(feature, qualifier_key, prefix):
    """
    Generic function to get the first BioPython feature qualifier value that starts with a given prefix.
    For example to get 'PF01234' from qualifiers of {'db_xref': ['PFAM: PF01234']}
    you can use get_qualifier_without_prefix(feature, 'db_xref', 'PFAM: ')
    :param feature: SeqIO feature to get qualifier value from
    :param qualifier_key: Key of qualifier to use
    :param prefix: Prefix to look for
    :return: Qualifier value after given prefix
    """
    qualifiers = [q for q in feature.qualifiers.get(qualifier_key) if prefix in q]
    return qualifiers[0].replace(prefix, '') if qualifiers else None

def get_domain_cf_probability(domain_feature):
    return float(get_qualifier_without_prefix(domain_feature, 'note', 'ClusterFinder probability: '))

def get_domain_accession(domain_feature):
    return get_qualifier_without_prefix(domain_feature, 'db_xref', 'PFAM: ')

def get_first_qualifier(feature, qualifier_key):
    """
    Get first value of given qualifier or None if it is not present.
    :param feature: Feature to get qualifier from
    :param qualifier_key: Qualifier key to get
    :return: first value of given qualifier or None if it is not present
    """
    return feature.qualifiers.get(qualifier_key)[0] if feature.qualifiers.get(qualifier_key) else None

def get_features_by_locus(record, feature_type):
    """
    Get dict of features of given type by their locus_tag qualifier
    :param record: Record to extract features from
    :param feature_type: Type of features to consider
    :return: dict of features of given type by their locus_tag qualifier
    """
    features = {}
    for feature in record.features:
        if feature.type == feature_type:
            locus_tag = get_first_qualifier(feature, 'locus_tag')
            if locus_tag in features:
                raise AttributeError('Locus tag "{}" found more than once in {}'.format(locus_tag, record.id))
            features[locus_tag] = feature
    return features

def antismash_cf_prediction(records):
    """
    Convert the final antiSMASH GenBank file with ClusterFinder predictions into a DataFrame of protein domains
    :param records: SeqIO records from the final antiSMASH GenBank file
    :return: DataFrame of protein domains with ClusterFinder 'prediction' column
    """
    domains = []
    for record in records:
        print('Processing record {}'.format(record.id))
        proteins_by_locus = get_features_by_locus(record, 'CDS')
        genes_by_locus = get_features_by_locus(record, 'gene')
        for feature in record.features:
            if feature.type == 'PFAM_domain':
                locus_tag = get_first_qualifier(feature, 'locus_tag')
                protein = proteins_by_locus.get(locus_tag)
                gene = genes_by_locus.get(locus_tag)
                gene_start = int(gene.location.start)
                gene_end = int(gene.location.end)
                gene_strand = int(gene.location.strand)
                protein_id = get_first_qualifier(protein, 'protein_id') or locus_tag
                if not protein_id:
                    raise AttributeError('No protein_id or locus tag for feature {}'.format(protein))
                domains.append({
                    'prediction': get_domain_cf_probability(feature),
                    'pfam_id': get_domain_accession(feature),
                    'evalue': float(get_first_qualifier(feature, 'evalue')),
                    'gene_start': gene_start,
                    'gene_end': gene_end,
                    'gene_strand': gene_strand,
                    'locus_tag': locus_tag,
                    'protein_id': protein_id,
                    'contig_id': record.id
                })
    return pd.DataFrame(domains)[['contig_id','gene_start','gene_end','gene_strand','locus_tag','protein_id','pfam_id','evalue','prediction']]


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Path to AntiSMASH final GenBank file.", metavar="FILE")
    parser.add_argument("-t", "--true", dest="true", required=False,
                        help="Path to true output CSV file.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    records = SeqIO.parse(options.input, 'genbank')
    domains = antismash_cf_prediction(records)

    if options.true:
        print('Adding in_cluster column with true output from {}'.format(options.true))
        true_output = pd.read_csv(options.true).groupby("protein_id")["in_cluster"].mean().to_dict()
        domains["in_cluster"] = domains["protein_id"].apply(lambda p: true_output.get(p)).fillna(method="ffill")

    domains.to_csv(options.output, index=False)

    print('Saved prediction of {} domains to {}'.format(len(domains), options.output))

