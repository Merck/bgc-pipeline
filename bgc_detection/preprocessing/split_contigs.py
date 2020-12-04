#!/usr/bin/env python
# David Prihoda
# Split one Genbank/FASTA/Embl file with multiple contigs into one file per contig

import argparse
import os
from Bio import SeqIO

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--format", dest="format", required=True,
                      help="Input file format (genbank/fasta/embl).", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output folder name.", metavar="FILE")
    parser.add_argument(dest="input", nargs='+',
                      help="Input genome file paths.", metavar="FILES")

    options = parser.parse_args()

    os.makedirs(options.output)

    num = 0
    for path in options.input:
        contigs = SeqIO.parse(path, options.format)
        _, extension = os.path.splitext(path)

        for contig in contigs:
            out_path = os.path.join(options.output, contig.id + extension)
            print("Saving {} to {}".format(contig.id, out_path))
            if os.path.exists(out_path):
                raise ValueError('Duplicate contig version, file already exists: '+out_path)
            SeqIO.write(contig, out_path, options.format)
            num += 1

    print('Saved {} contigs to {}.'.format(num, options.output))