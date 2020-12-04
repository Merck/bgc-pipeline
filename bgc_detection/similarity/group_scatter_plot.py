#!/usr/bin/env python
# David Prihoda
# Create a scatterplot from a CSV file with numeric columns
# Can visualize multi-label classes using point color and outline

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import manifold
import itertools
import seaborn as sns
import numpy as np

COLORS = [
    'blue',
    'deeppink',
    'green',
    'darkcyan',
    'indigo',
    'orangered',
    'brown',
    'gold',
    'darkkhaki',
    'lime',
    'lightcoral',
    'turquoise'
]

REMAINING_LABEL = '- remaining -'
MISSING_LABEL = '- missing -'

PALETTES = {
    'PRODUCT_CLASS': {
        'Polyketide' : 'deeppink',
        'NRP' : 'lime',
        'RiPP' : 'blue',
        'Other' : 'darkcyan',
        'Saccharide' : 'orange',
        'Terpene': 'turquoise',
        'Alkaloid' : 'purple',
        '?' : None,
        'Nucleoside': 'brown',
        'Geneswap_Negative' : 'black',
        'Candidate' : 'grey',
        'MIBiG_Positive' : '#2255ee',
        REMAINING_LABEL : 'grey',
        MISSING_LABEL : 'black'
    },
    'CHEM_ACT' : {
        'Antibacterial': 'deeppink',
        'Antifungal': 'orange',
        'Cytotoxic': 'blue',
        'Inhibitor': 'darkcyan',
        'Other': 'purple',
        'Signalling': 'turquoise',
        'Surfactant': 'lime',
        'Unknown': 'brown',
        REMAINING_LABEL : 'grey',
        MISSING_LABEL : 'black'
    }
}
def group_scatter_plot(samples, figure_size=10, groupby=None, legend_loc="lower left", verbose=1, linewidth=2, legend_hybrids=True, axes_label=None, swap_axes=False, **kwargs):
    """
    Create scatterplot from the first two numeric values in given DataFrame, grouped and colored by given column
    :param samples: DataFrame with samples, each sample has numeric columns plus optionally a groupby column and 'color', 'size', 'marker' and 'alpha'.
    :param figure_size: Width=Height of produced figure
    :param groupby: Column to use for grouping and coloring
    :param legend_loc: Legend location
    :param verbose: Verbosity
    :param linewidth: Point Line width
    :param legend_hybrids: Whether to show legend for multi-valued groups
    :param axes_label: Label of each axis, will produce 'label1' for x-axis and 'label2' for y-axis
    :param swap_axes: Swap x and y axis
    :param kwargs: Additional arguments to pass to scatterplot function
    :return: Figure with scatterplot
    """
    samples = samples.copy()
    fig, ax = plt.subplots(1, 1, figsize=(figure_size, figure_size))
    is_grouped = True
    if not groupby:
        is_grouped = False
        samples[groupby] = None
    if 'color' not in samples:
        samples['color'] = '#1144ee'
    if 'size' not in samples:
        samples['size'] = 10
    if 'marker' not in samples:
        samples['marker'] = 'o'
    if 'alpha' not in samples:
        samples['alpha'] = 0.7

    numeric_cols = list(samples.select_dtypes(include=[np.number]).columns)
    print('Numeric: ', numeric_cols)
    axis_index = (0, 1) if not swap_axes else (1, 0)
    for group_args, group_samples in samples.groupby(['marker', groupby, 'color', 'alpha']):
        marker, group, color, alpha = group_args
        if is_grouped and not group:
            print('Skipping {} samples of blank group: {}'.format(len(group_samples), group_args))
            continue
        edgecolors = 'none'
        label = group
        if ',' in color:
            color, edgecolors = color.split(',')
            # Remove marker for edges
            if ':' in edgecolors:
                edgecolors = edgecolors.split(':')[0]
            if not legend_hybrids:
                label = ''
        if ':' in color:
            color, marker = color.split(':')

        if verbose:
            print('Drawing {}x {}'.format(len(group_samples), group))
            print(' ', group_args)
            print(group_samples.head())
            print('-' * 80)
        group_samples.plot(kind='scatter', x=numeric_cols[axis_index[0]], y=numeric_cols[axis_index[1]], label=label, ax=ax, edgecolors=edgecolors,
                               s=group_samples['size'], marker=marker, linewidth=linewidth, alpha=alpha, color=color, **kwargs)

    if axes_label:
        ax.set_xlabel('{}{}'.format(axes_label, axis_index[0] + 1))
        ax.set_ylabel('{}{}'.format(axes_label, axis_index[1] + 1))

    bbox = None
    if 'outside' in legend_loc:
        legend_loc = legend_loc.replace('outside', '').strip()
        bbox = (1, 1)
        if 'lower' in legend_loc:
            bbox = (1, 0)
    lgnd = plt.legend(loc=legend_loc, bbox_to_anchor=bbox, scatterpoints=1, fontsize=10, frameon=False)
    if lgnd:
        if not legend_loc:
            lgnd.remove()
        else:
            for handle in lgnd.legendHandles:
                handle._sizes = [40]
                handle.set_alpha(1)
    return fig

def default(l, idx, default):
  try:
    return l[idx]
  except IndexError:
    return default
  
def get_top_groups(groups, num_groups):
    counts = groups.value_counts()
    if not num_groups:
        return counts.index.values
    print('Groups:')
    print(counts)
    return counts.index.values[:num_groups]

def apply_group_color(group, palette, default='grey'):
    if isinstance(group, str) and ';' in group:
        first_group, second_groups = group.split(';', 1)
        return palette.get(first_group, default) + ',' + palette.get(second_groups, default)
    return palette.get(group, default)

def add_group_properties(samples, groups, palette=None, num_groups=None, color_offset=0, **kwargs):
    top_groups = get_top_groups(groups, num_groups=num_groups)
    groups = groups.apply(lambda g: g if g in top_groups else (MISSING_LABEL if g is np.nan or not g else REMAINING_LABEL))
    if not palette:
        palette = {g: COLORS[(i+color_offset) % len(COLORS)] for i, g in enumerate(top_groups)}
        palette[MISSING_LABEL] = 'black'
        palette[REMAINING_LABEL] = 'grey'
    # Assign groups and properties to given samples
    groups = groups.loc[np.intersect1d(groups.index, samples.index)]
    samples.loc[groups.index, 'group'] = groups
    samples.loc[groups.index, 'color'] = groups.apply(lambda g: apply_group_color(g, palette))
    for prop, val in kwargs.items():
        if val is not None:
            print('Adding {}: {} = {}'.format(len(groups), prop, val))
            samples.loc[groups.index, prop] = val
    return top_groups

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="Reduced feature matrix CSV file.", metavar="FILE")

    parser.add_argument("--group", dest="groups", action='append', default=[],
                      help="Group properties CSV file.", metavar="FILE")

    parser.add_argument("--index", dest="index", action='append', default=[],
                      help="Index column in properties table.", metavar="FILE")

    parser.add_argument("--label", dest="label", default=[], action='append',
                      help="Subgroup column in properties table if present, otherwise group label.", metavar="STRING")

    parser.add_argument("--label-prefix", dest="label_prefix", action='append', default=[],
                      help="Add a prefix to label value.", metavar="STR")

    parser.add_argument("--marker", dest="marker", default=[], action='append',
                      help="Group marker.", metavar="STRING")

    parser.add_argument("--alpha", dest="alpha", default=[], type=float, action='append',
                      help="Marker color alpha (opacity).", metavar="FLOAT")

    parser.add_argument("--size", dest="size", default=[], type=float, action='append',
                      help="Marker size.")

    parser.add_argument("--min-size", dest="min_size", default=[], type=float, action='append',
                      help="Minimum marker size to add to size column.")

    parser.add_argument("--size-column", dest="size_column", default=[], action='append',
                      help="Marker size multiplication column.")

    parser.add_argument("--num-groups", dest="num_groups", default=[], type=int, action='append',
                      help="Number of top groups to color.", metavar="INT")

    parser.add_argument("--figure-size", dest="figure_size", default=10, type=float,
                      help="Figure size.", metavar="FLOAT")

    parser.add_argument("--xlim-from", dest="xlim_from", default=None, type=float,
                      help="X axis limit from.", metavar="FLOAT")

    parser.add_argument("--xlim-to", dest="xlim_to", default=None, type=float,
                      help="X axis limit to.", metavar="FLOAT")

    parser.add_argument("--palette", dest="palette",
                      help="Color palette: ({}).".format(PALETTES.keys()), metavar="STRING")

    parser.add_argument("--legend-no-hybrids", dest="legend_hybrids", action="store_false", default=True,
                      help="Hide hybrids (multiple classes) from legend.")

    parser.add_argument("--legend", dest="legend", default='lower left',
                      help="Legend location ([outside] lower/upper left/right).", metavar="STRING")

    parser.add_argument("--title", dest="title",
                      help="Plot title.", metavar="STRING")

    parser.add_argument("--ax", dest="ax",
                      help="Axis labels.", metavar="STRING")

    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output png file path.", metavar="FILE")

    parser.add_argument("--swap-axes", dest="swap_axes", action="store_true",
                      help="Swap axes of tsne plot.")

    parser.add_argument("-c", "--csv", dest="csv", required=False,
                      help="Output csv file path.", metavar="FILE")

    options = parser.parse_args()

    samples = pd.read_csv(options.input)
    if options.groups:
        samples.set_index(samples.columns[0], inplace=True)
        # Read all property files indexed by configured column
        groups = [pd.read_csv(path).set_index(default(options.index, i, 'contig_id')) for i, path in enumerate(options.groups)]
    else:
        # Handle all samples as one single group
        groups = [samples[[]]]
    samples['group'] = None

    palette = PALETTES[options.palette] if options.palette else None
    min_size = 0
    size = 10
    alpha = 0.7
    marker = 'o'
    num_groups = None
    color_offset = 0
    for i, group_properties in enumerate(groups):
        if len(options.min_size) > i:
            min_size = options.min_size[i]
        if len(options.size) > i:
            size = options.size[i]
        if len(options.alpha) > i:
            alpha = options.alpha[i]
        if len(options.marker) > i:
            marker = options.marker[i]
        if len(options.num_groups) > i:
            num_groups = options.num_groups[i]
        # Get parameters or defaults
        label = default(options.label, i, 'Vector')
        size_column = default(options.size_column, i, None)
        size_multiplier = group_properties[size_column] if size_column else 1
        sizes = size * size_multiplier + min_size if size else None
        if label in group_properties:
            subgroups = group_properties[label]#.astype(np.str)
            prefix = default(options.label_prefix, i, "")
            subgroups = subgroups.apply(lambda label: prefix+label.replace(';', ';'+prefix) if isinstance(label, str) else label)
        else:
            subgroups = pd.Series(label, group_properties.index)
        added_groups = add_group_properties(
            samples=samples,
            groups=subgroups,
            palette=palette,
            num_groups=num_groups,
            color_offset=color_offset,
            size=sizes,
            alpha=alpha,
            marker=marker
        )
        color_offset += len(added_groups)

    fig = group_scatter_plot(
        samples=samples,
        groupby='group',
        figure_size=options.figure_size,
        legend_loc=options.legend,
        title=options.title,
        legend_hybrids=options.legend_hybrids,
        axes_label=options.ax,
        swap_axes=options.swap_axes
    )
    if options.xlim_from or options.xlim_to:
        fig.axes[0].set_xlim(options.xlim_from, options.xlim_to)
    fig.savefig(options.output, dpi=200, bbox_inches='tight')
    print('Saved plot to:', options.output)
    if options.csv:
        samples.to_csv(options.csv, index=False)
        print('Saved csv to:', options.csv)


