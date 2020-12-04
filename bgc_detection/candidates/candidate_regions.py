#!/usr/bin/env python
# David Prihoda
# Plot the "region plot" of BGC candidates in a bacterial genomes (horizontal colored lines for each model).

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def candidate_regions(cands, safety_limit=50, xlim=0, xstep=100000, colors=None):
    """
    Plot the "region plot" of BGC candidates in a bacterial genomes (horizontal colored lines for each model).
    :param cands: DataFrame of candidates with 'contig_id' column defining the contig and 'model' defining the model used to predict the candidate
    :param safety_limit: Maxmimum number of contigs to plot
    :param xlim: Starting genomic coordinate
    :param xstep: Step for X-axis ticks
    :param colors: Dict of colors by 'model'
    :return: Figure for the "region plot" of BGC candidates
    """
    cands_by_contig = cands.groupby('contig_id')
    if len(cands_by_contig) > safety_limit:
        raise AttributeError('You probably did not want to plot more than {} contigs! '
                             'Otherwise increase the safety_limit parameter.'.format(safety_limit))

    models = cands['model'].unique()

    lengths = cands_by_contig['nucl_end'].max()
    max_end = lengths.max()

    with plt.style.context(('ggplot')):
        rows = len(cands_by_contig)
        num_models = len(models)
        width = 15 + (max_end-xlim) / 230000
        fig, axes = plt.subplots(rows, 1, figsize=(width, 2 + 0.20 * rows * num_models))
        if rows == 1:
            axes = [axes]
        for i, (contig_id, contig_cands) in enumerate(cands_by_contig):
            print('{} of length {}'.format(contig_id, lengths[contig_id]))
            ax = axes[i]
            ax.set_title(contig_id, loc='left')
            ax.set_facecolor('white')
            ax.set_yticks([])
            ax.set_xlim([xlim, max_end])
            ax.set_ylim([-0.1, num_models + 0.5])
            ax.set_xticks(range(xlim, int(lengths[contig_id]), xstep))
            ax.set_xticklabels(
                ['{:.0f}kb'.format(x / 1e3) if x < 1e6 else '{:.1f}Mb'.format(x / 1e6) for x in ax.get_xticks()])
            cands_by_model = contig_cands.groupby('model')
            for m, model in enumerate(models):
                if model not in cands_by_model.groups:
                    continue
                model_cands = cands_by_model.get_group(model)

                if model == 'true_output':
                    color = 'grey'
                    for c, cand in model_cands.iterrows():
                        ax.axvspan(cand['nucl_start'], cand['nucl_end'], color='black', alpha=0.13)
                else:
                    color = colors[model] if colors else plt.cm.tab10(m)

                x = model_cands[['nucl_start', 'nucl_end']].values.reshape(1, -1)[0]
                y = np.ones(x.shape) * (num_models - m)
                y[1::2] = np.nan
                for d in np.arange(-0.08, 0.08, 0.01):
                    ax.step(x, y + d, color=color, where='post', lw=0.25, label=None)
                # hidden plot for legend
                ax.step([-100, -120], [0, 0], color=color, where='post', lw=4, label=model)

        axes[0].legend(loc='upper right', bbox_to_anchor=(-0.01, 1.0))
        fig.patch.set_facecolor('white')
        fig.tight_layout()
        fig.subplots_adjust(left=0.03 + 3 / width)

    return fig


def create_and_save(cands, path, **kwargs):
    fig = candidate_regions(cands, **kwargs)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print('Saved regions plot to: ', path)


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--true", dest="true", required=False,
                      help="True BGCs candidate csv file path.", metavar="FILE")
    parser.add_argument("--candidates", dest="candidates", required=True, action='append',
                      help="Model candidates CSV file path.", metavar="FILE")
    parser.add_argument("--color", dest="colors", action='append', required=True,
                      help="Model color.", metavar="STRING")
    parser.add_argument("--label", dest="labels", action='append', required=True,
                      help="Model label.", metavar="STRING")
    parser.add_argument("-s", "--separate", dest="separate", action='store_true', default=False,
                        help="Create one separate file for each contig.")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output file path.", metavar="FILE")

    options = parser.parse_args()

    cands = []

    if options.true:
        true_output = pd.read_csv(options.true)
        true_output['model'] = 'true_output'
        cands.append(true_output)

    colors = {}
    for i, path in enumerate(options.candidates):
        name = options.labels[i]
        c = pd.read_csv(path)
        c['model'] = name
        if options.colors and len(options.colors) > i:
            colors[name] = options.colors[i]
        cands.append(c)
    cands: pd.DataFrame = pd.concat(cands)

    if options.separate:
        print('Creating output directory {}'.format(options.output))
        os.mkdir(options.output)

        grouped = cands.groupby('contig_id')
        print('Drawing {} contigs...'.format(len(grouped)))
        for contig_id, contig_cands in grouped:
            path = os.path.join(options.output, '{}.png'.format(contig_id))
            create_and_save(contig_cands, path, colors=colors)
    else:
        print('Drawing region plot...')
        create_and_save(cands, options.output, colors=colors)

