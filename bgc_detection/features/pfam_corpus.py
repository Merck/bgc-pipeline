#!/usr/bin/env python
# David Prihoda
# Generate pfam corpus for training pfam2vec from Domain CSV files
# Will produce a single line for each input file, with pfam IDs separated by spaces

import argparse

import os
import pandas as pd
import glob

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Input folder with domain csv files.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output corpus txt file.", metavar="FILE")
    parser.add_argument("-e", "--maxevalue", dest="maxevalue", required=True, type=float,
                        help="Maximum domain independent e-value.", metavar="FLOAT")
    options = parser.parse_args()

    count = 0
    paths = list(glob.glob(os.path.join(options.input, '*.csv')))
    for path in paths:
        count += 1
        with open(options.output, 'a') as corpusfile:
            print("{} ({}/{})".format(path, count, len(paths)))
            domains = pd.read_csv(path)
            domains = domains[domains['evalue'] <= options.maxevalue]
            corpusfile.write(' '.join(domains['pfam_id'].values))
            corpusfile.write('\n')

    print('Saved corpus of {} domain files to {}'.format(count, options.output))

