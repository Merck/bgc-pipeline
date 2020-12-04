#!/usr/bin/env python
# David Prihoda
# Extract CSV with BGC species and biosynthetic class information from GenBank and JSON MIBiG files

import argparse
import pandas as pd
from Bio import SeqIO
import json
import os
import glob

def get_cluster_ids(gbk_folder):
    return [os.path.basename(file.replace('.gbk', '')) for file in glob.glob(os.path.join(gbk_folder,'*.gbk'))]

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
    print('Getting properties for {} clusters...'.format(len(cluster_ids)))
    properties = []
    for cluster_id in cluster_ids:
        try:
            with open(os.path.join(options.json, cluster_id+'.json'), 'r') as jsonfile:
                json_props = json.load(jsonfile)
                classes = json_props.get('general_params').get('biosyn_class')
                c_class = ';'.join(classes)
        except:
            print('No json properties for cluster {}, using "?"'.format(cluster_id))
            c_class = '?'
        records = list(SeqIO.parse(os.path.join(options.gbk, cluster_id+'.gbk'), 'genbank'))
        for record in records:
            source = record.annotations.get('source')
            proteins = [f for f in record.features if f.type == 'CDS']
            properties.append({
                'BGC_ID' : cluster_id,
                'contig_id' : record.id,
                'num_proteins': len(proteins),
                'classes' : c_class,
                'source' : source,
                'species' : source.split(' ')[0]
            })
    properties = pd.DataFrame(properties)
    properties[['BGC_ID','contig_id','num_proteins','classes','source','species']].to_csv(options.output, index=False)

    print('Saved to', options.output)

