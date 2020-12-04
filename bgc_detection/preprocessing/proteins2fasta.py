#!/usr/bin/env python
# David Prihoda
# Convert proteins from an annotated GenBank file into a FASTA file

import argparse
import glob

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Input genome file/folder path.", metavar="FILE")
    parser.add_argument("-e", "--extension", dest="extension", default=None,
                      help="Include all files matching extension.", metavar="text")
    parser.add_argument("-f", "--format", dest="format", required=True,
                      help="Input file format (genbank/embl).", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output fasta file name.", metavar="FILE")
    options = parser.parse_args()

    with open(options.output, 'a') as out:
        print('Appending to', options.output)
        inputs = [options.input]
        if options.extension:
            inputs = glob.glob(os.path.join(options.input, '*'+options.extension))
        for input in inputs:
            print('Adding file {}'.format(input))
            for record in SeqIO.parse(input, options.format):
                print('Converting proteins from {}...'.format(record.id))
                proteins = []
                for feature in record.features:
                    if feature.type == 'CDS' and feature.qualifiers.get('translation'):
                        locus_tags = feature.qualifiers.get('locus_tag')
                        locus_tag = str(locus_tags[0]) if locus_tags else 'unknown_locus_tag'
                        # a unique protein id is required, if protein_id is not present, we need a locus tag
                        location = str(feature.location.start) + '-' + str(feature.location.end)
                        protein_ids = feature.qualifiers.get('protein_id')
                        if protein_ids:
                            protein_id = str(protein_ids[0])
                        elif locus_tags:
                            protein_id = str(locus_tags[0])
                        else:
                            protein_id = 'unknown_protein'+location
                        r_id = record.id + '|' + locus_tag + '|' + protein_id + '|' + location + '|' + str(feature.location.strand)
                        translation = Seq(feature.qualifiers.get('translation')[0], generic_protein)
                        r = SeqRecord(translation, r_id, description='')
                        proteins.append(r)

                SeqIO.write(proteins, out, 'fasta')

        print('Saved to', options.output)

