#!/usr/bin/env python
# David Prihoda
# Predict BGC scores for a Domain CSV file using a trained model

from utils import io
from pipeline import PipelineWrapper
import argparse
import os
from multiprocessing import Pool
import threading
import pandas as pd

from candidates.average_protein_prediction import average_protein_prediction

cache = threading.local()

def run_prediction(domains, pipeline, whole=False):
    """
    Get BGC prediction score for given Domain DataFrame, add it as 'prediction' column.
    :param domains: Domain DataFrame, multiple samples marked by different 'contig_id' will be predicted separately
    :param pipeline: Trained BGC detection Pipeline
    :param whole: Always predict whole sequence together, even if different 'contig_ids' are present.
    :return: Original Domain DataFrame with 'prediction' column added.
    """
    # Single sample
    if whole or 'contig_id' not in domains.columns:
        print('Predicting single sample (contig_id column not present).')
        domains['prediction'] = pipeline.predict(domains)
        return domains
    # Multiple samples
    else:
        samples = io.domains_to_samples(domains, 'contig_id')
        print('Predicting {} samples...'.format(len(samples)))
        predictions = []
        for sample in samples:
            prediction = sample.copy()
            prediction['prediction'] = pipeline.predict(sample)
            predictions.append(prediction)

        merged: pd.DataFrame = pd.concat(predictions)
        return merged

def run_prediction_task(args):
    i, path, options = args
    pipeline = getattr(cache, 'pipeline', None)
    if not pipeline:
        print('Creating model', MODEL_PATH)
        pipeline = PipelineWrapper.load(MODEL_PATH)
        cache.pipeline = pipeline

    domains = io.read_domains(path, options.maxevalue, None)
    prediction = run_prediction(domains, pipeline, whole=options.whole)
    if options.avg:
        prediction = average_protein_prediction(prediction, prediction['prediction'])
    output_path = options.output
    if not output_path:
        filename = os.path.basename(path)
        output_path = os.path.join(options.dir, filename)
    print('Saving prediction #{} to {}'.format(i + 1, output_path))
    prediction.to_csv(output_path, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--model", dest="model", required=True,
                        help="Path to trained model pickle file.", metavar="FILE")
    parser.add_argument("-d", "--dir", dest="dir", required=False,
                        help="Output dir path.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=False,
                        help="Output file path.", metavar="FILE")
    parser.add_argument("-e", "--maxevalue", dest="maxevalue", required=True, type=float,
                        help="Maximum domain independent e-value.", metavar="FLOAT")
    parser.add_argument("-a", "--avg", dest="avg", action='store_true',
                        help="Average protein prediction.")
    parser.add_argument("--whole", dest="whole", action='store_true',
                        help="Discard contig_id information and predict whole sequence at once.")
    parser.add_argument(dest='samples', nargs='+',
                        help="Paths to samples to predict.", metavar="SAMPLES")
    options = parser.parse_args()

    if (options.dir and options.output) or (not options.dir and not options.output):
        raise AttributeError('Specify output directory using --dir or output file using --output.')

    MODEL_PATH = options.model

    tasks = [(i, path, options) for i, path in enumerate(options.samples)]

    pool = Pool()
    pool.map(run_prediction_task, tasks)
    print('Done.')


