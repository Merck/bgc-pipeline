#!/usr/bin/env python
# David Prihoda
# Create pfam2vec embedding from corpus of pfam IDs
# Corpus words should be separated by spaces, documents (genomes) should be separated by newlines

import word2vec
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True, action='append',
                      help="Input corpus txt file path (repeat for multiple files). "
                           "Words should be separated by spaces, documents (genomes) should be separated by newlines.", metavar="FILE")
    parser.add_argument("-d", "--dimensions", dest="dimensions", required=True, action='append', type=int,
                      help="Number of dimensions.", metavar="FILE")
    parser.add_argument("-w", "--window", dest="window", type=int, default=5,
                      help="Window size.", metavar="FILE")
    parser.add_argument("-it", "--iter", dest="iter", type=int, action='append',
                      help="Number of iterations.", metavar="FILE")
    parser.add_argument("-m", "--method", dest="method", action='append',
                      help="Method (cbow, skipgram). Repeat for more methods.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output folder.", metavar="FILE")
    options = parser.parse_args()

    METHOD_CODES = {
        'skipgram': 0,
        'cbow': 1
    }

    if not options.iter:
        options.iter = [5]

    if not options.method:
        options.method = ['skipgram', 'cbow']

    for corpus_path in options.input:
        for num_features in options.dimensions:
            for iters in options.iter:
                for method in options.method:
                    corpus_name = os.path.splitext(os.path.basename(corpus_path))[0]
                    out_path = os.path.join(options.output, 'pfam2vec_{}_{}_{}dim_{}win_{}iter.bin'.format(corpus_name, method, num_features, options.window, iters))
                    print('-'*80)
                    print('Building Word2Vec:{} model for {}, producing {} dimensions'.format(method, corpus_name, num_features))
                    word2vec.word2vec(corpus_path, out_path, size=num_features, cbow=METHOD_CODES[method], iter_=iters, window=options.window, verbose=True)
                    print('\nSaved to ', out_path)
