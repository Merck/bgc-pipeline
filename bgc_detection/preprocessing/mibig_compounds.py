#!/usr/bin/env python
# David Prihoda
# Extract compound information CSV from GenBank and JSON MIBiG files

import argparse
import pandas as pd
from Bio import SeqIO
import json
import os
from mibig_properties import get_cluster_ids

COMPOUND_ACTIVITY = ['Antibacterial', 'Antifungal', 'Antioxidant', 'Antiviral',
                     'Cytotoxic', 'Inhibitor', 'Other', 'Pigment', 'Siderophore', 'Signalling', 'Surfactant', 'Unknown']


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gbk", dest="gbk", required=True,
                      help="Input cluster gbk folder.", metavar="FILE")
    parser.add_argument("-j", "--json", dest="json", required=True,
                      help="Input cluster json folder.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output csv file name.", metavar="FILE")
    options = parser.parse_args()

    cluster_ids = get_cluster_ids(options.gbk)
    print('Getting compound properties for {} clusters...'.format(len(cluster_ids)))
    properties = []
    for cluster_id in cluster_ids:
        try:
            with open(os.path.join(options.json, cluster_id+'.json'), 'r') as jsonfile:
                json_props = json.load(jsonfile)
        except:
            print('No json properties for cluster {}'.format(cluster_id))
            continue

        compounds = json_props.get('general_params').get('compounds', [])
        compound_names = [c.get('compound') for c in compounds]
        activities = sorted(list(set([a for c in compounds for a in c.get('chem_act', []) if a != ''])))
        props = {
            'BGC_ID': cluster_id,
            'num_compounds': len(compound_names),
            'compounds': ';'.join(compound_names),
            'num_activity': len(activities),
            'activity': ';'.join(activities),
            'any_known_activity': any([a != 'Unknown' for a in activities])
        }
        for chem_act in COMPOUND_ACTIVITY:
            props['any_' + chem_act.lower()] = int(any([chem_act in c.get('chem_act', []) for c in compounds]))

        records = list(SeqIO.parse(os.path.join(options.gbk, cluster_id + '.gbk'), 'genbank'))
        for record in records:
            record_props = dict(props)
            record_props['contig_id'] = record.id
            properties.append(record_props)

    columns = ['BGC_ID', 'contig_id', 'num_compounds', 'compounds', 'num_activity', 'activity', 'any_known_activity']
    columns += ['any_'+act.lower() for act in COMPOUND_ACTIVITY]
    properties = pd.DataFrame(properties)[columns]
    print(properties.head())
    print('Total:')
    print(properties.sum())

    properties.to_csv(options.output, index=False)
    print('Saved to', options.output)

