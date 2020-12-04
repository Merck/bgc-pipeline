#!/usr/bin/env python
# David Prihoda
# Convert a multi-valued CSV column into binary columns for each unique value
# Used to prepare multi-label input based on 'classes' column with zero to multiple classes separated by ';'.

import argparse
import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer

def classes_to_multilabel(properties, index, column, delimiter=';'):
    """
    Create binary multi-label DataFrame with one column for each unique value
    in given column of given DataFrame, indexed by a second given column.
    The column will be equal to 0 if the value was not present and 1 if it was.
    :param properties: DataFrame with the multi-valued column and the index column
    :param index: Name of column to use as index
    :param column: Name of multi-valued column to get values from
    :param delimiter: Delimiter of the multi-valued column
    :return: Binary multi-label DataFrame with one column for each unique value
    in given column, indexed by a second given column.
    """
    class_sublists = properties[column].apply(lambda c: c.split(delimiter))
    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(class_sublists)
    ml = pd.DataFrame(Y, index=properties[index], columns=mlb.classes_)
    print('Properties:')
    print(properties.head())
    print('Multi-label:')
    print(ml.head())
    return ml


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="Properties csv file.", metavar="FILE")
    parser.add_argument("--index", dest="index", required=False, default='contig_id',
                        help="Column to use as index.", metavar="STRING")
    parser.add_argument("--column", dest="column", required=True,
                        help="Column to split as multilabel output.", metavar="STRING")
    parser.add_argument("--delimiter", dest="delimiter", required=False, default=';',
                        help="Delimiter to use for splitting the column value.", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    properties = pd.read_csv(options.input)

    vectors = classes_to_multilabel(properties, index=options.index, column=options.column, delimiter=options.delimiter)

    print('Saving {} multilabel outputs to: {}'.format(len(vectors), options.output))
    vectors.to_csv(options.output, index=True)
    print('Done.')

