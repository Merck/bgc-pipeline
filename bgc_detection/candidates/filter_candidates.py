#!/usr/bin/env python
# David Prihoda
# Filter a BGC candidate CSV file by different fields - coverage, nucleotide length or number of domains
# Used e.g. to get novel candidates (with 0 coverage)

import argparse
import pandas as pd

def get_top_candidates(target_cands, max_coverage=None, min_nucleotides=None, min_domains=None):
    """
    Filter BGC candidate DataFrame by different fields - coverage, nucleotide length or number of domains
    :param target_cands: DataFrame of BGC candidates
    :param max_coverage: maximum allowed coverage
    :param min_nucleotides: minimum allowed number of nucleotides
    :param min_domains: minimum number of domains
    :return: BGC candidate DataFrame filtered by given fields
    """
    top_cands = target_cands.sort_values(by='nucl_length', ascending=False)
    if max_coverage is not None:
        prev_size = len(top_cands)
        top_cands = top_cands[top_cands['coverage'] <= max_coverage]
        print('Filtered candidates by maximum coverage of {:2f}, {}/{} remained.'.format(max_coverage, len(top_cands), prev_size))
    if min_nucleotides is not None:
        prev_size = len(top_cands)
        top_cands = top_cands[top_cands['nucl_length'] >= min_nucleotides]
        print('Filtered candidates by minimum of {} nucleotides, {}/{} remained.'.format(min_nucleotides, len(top_cands), prev_size))
    if min_domains is not None:
        prev_size = len(top_cands)
        top_cands = top_cands[top_cands['num_all_domains'] >= min_domains]
        print('Filtered candidates by minimum of {} domains, {}/{} remained.'.format(min_domains, len(top_cands), prev_size))
    return top_cands

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Target model candidate csv file path.", metavar="FILE")
    parser.add_argument("-c", "--max-coverage", dest="max_coverage", required=False, type=float,
                      help="Maximum allowed coverage by other models.", metavar="FILE")
    parser.add_argument("-n", "--min-nucleotides", dest="min_nucleotides", required=False, type=int,
                      help="Minimum required nucleotide length.", metavar="FILE")
    parser.add_argument("-d", "--min-domains", dest="min_domains", required=False, type=int,
                      help="Minimum required pfam domains.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    target_cands = pd.read_csv(options.input)

    top_cands = get_top_candidates(
        target_cands=target_cands,
        max_coverage=options.max_coverage,
        min_nucleotides=options.min_nucleotides,
        min_domains=options.min_domains
    )

    top_cands.to_csv(options.output, index=False)
    print('Saved {} top candidates to: {}'.format(len(top_cands), options.output))

