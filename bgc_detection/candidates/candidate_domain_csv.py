#!/usr/bin/env python
# David Prihoda
# Create a minimal Domain CSV file from the Candidate CSV file using the 'pfam_ids' column

import argparse
import pandas as pd

def split_on_column(df, index_column, split_column, new_index_column, new_split_column, split_char=';'):
    """
    Split a given DataFrame column (split_column) into multiple rows together with a second index column (index_column).
    :param df: DataFrame to split
    :param index_column: Column to include together with the split_column
    :param split_column: Column to split into multiple rows
    :param new_index_column: Name for index column
    :param new_split_column: Name for split column
    :param split_char: Character to split the column on
    :return: DataFrame where each value of split_column was split into multiple rows, where each row contains the same index_column value
    """
    # Produce a Series with index from index_column and values from split_column (split by split_char).
    # Each series constructor is called with (value, [split1, split2, split3]) which produces [(value, split1), (value, split2), (value, split3)]
    split = pd.concat([pd.Series(row[index_column], row[split_column].split(split_char))
               for _, row in df.iterrows()]).reset_index()
    # Add column names
    split.columns = [new_split_column, new_index_column]
    # Return reversed columns
    return split[[new_index_column, new_split_column]]

def candidate_domain_csv(cands, index_column):
    """
    Create a minimal Domain DataFrame from the Candidate CSV file using the 'pfam_ids' column
    :param cands: DataFrame of BGC candidates
    :param index_column: column to use as 'contig_id', should be 'candidate_hash'
    :return: a minimal Domain DataFrame from the Candidate CSV file using the 'pfam_ids' column
    """
    return split_on_column(cands, index_column, 'pfam_ids', 'contig_id', 'pfam_id')

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Model candidates CSV file path.", metavar="FILE")
    parser.add_argument("--index", dest="index", required=False, default='candidate_hash',
                      help="Index column.", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output domain CSV file path.", metavar="FILE")

    options = parser.parse_args()

    cands = pd.read_csv(options.input)
    domains = candidate_domain_csv(cands, index_column=options.index)

    domains.to_csv(options.output, index=False)
    print('Saved {} domains to: {}'.format(len(domains), options.output))

