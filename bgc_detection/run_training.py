#!/usr/bin/env python
# David Prihoda
# Train and save a BGC detection model using a JSON model config and a set of positive and negative set of samples - Domain CSVs.

try:
    from utils import io
    from pipeline import PipelineWrapper
except ModuleNotFoundError:
    from bgc_detection.utils import io
    from bgc_detection.pipeline import PipelineWrapper
import argparse
import time
import re
import os
import json

def sanitize_path(path):
    return re.sub('[^a-zA-Z0-9\-.]+', '_', path)

def create_progress_log_path(log_dir, model_output_path):
    date = int(time.time())
    label = sanitize_path(os.path.basename(model_output_path))
    filename = "{}-{}".format(date, label)
    return os.path.join(log_dir, filename)

def read_samples(sample_paths, evalue):
    all_samples = []
    all_y = []
    for sample_path in sample_paths or []:
        domains = io.read_domains(sample_path, max_evalue=evalue)
        samples, y_list = io.domains_to_samples(domains, 'contig_id', 'in_cluster')
        print('Loaded {} samples and {} domains from {}'.format(len(samples), len(domains), sample_path))
        all_samples += samples
        all_y += y_list
    return all_samples, all_y


def run_training(config, output_path, sample_paths, validation_sample_paths=None, evalue=None, progress_log_path=None, files=None, verbose=1):
    """
    Train a and save a BGC detection model using a JSON model config and a set of positive and negative set of samples - Domain DataFrames.
    :param config: Model config parsed from JSON
    :param output_path: Path where to save trained model pickle file
    :param sample_paths: List of paths to Domain CSV training files. Each Domain CSV can contain multiple sample sequences marked with the 'contig_id' column. Output is provided in 'in_cluster' column.
    :param validation_sample_paths: List of paths to validation samples.
    :param evalue: Use given evalue threshold, only domains with a lower evalue are considered.
    :param progress_log_path: Path to folder where to store logging files (e.g. TensorBoard).
    :param files: Dictionary of file paths to inject into the model config. For example "{myFolder}/dep.txt" can be replaced to "../my/value/dep.txt" using {"myFolder": "../my/value"} dictionary
    :param verbose: Verbosity
    """
    if files:
        pairs = files.items() if isinstance(files, dict) else files
        for key, path in pairs:
            config = replace_in_dict(config, "{"+key+"}", path)

    pipeline = PipelineWrapper.from_config(config)

    print('Reading train samples:')
    train_samples, train_y = read_samples(sample_paths, evalue=evalue)
    print('Reading validation samples:')
    validation_samples, validation_y = read_samples(validation_sample_paths, evalue=evalue)

    print('Progress will be saved to:', progress_log_path)
    pipeline.fit(
        samples=train_samples,
        y=train_y,
        debug_progress_path=progress_log_path,
        validation_samples=validation_samples,
        validation_y=validation_y,
        verbose=verbose
    )

    pipeline.save(output_path)

    print('-'*80)
    print('Trained model saved to:', output_path)
    print('Progress saved to:', progress_log_path)
    print('-'*80)

def replace_in_dict(d, key, value):
    if isinstance(d, dict):
        return {k: replace_in_dict(v, key, value) for k, v in d.items()}
    elif isinstance(d, list):
        return [replace_in_dict(v, key, value) for v in d]
    elif isinstance(d, str):
        return d.replace(key, value)
    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config", dest="config", required=True,
                        help="Path to JSON model config file.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Where to write trained model file.", metavar="FILE")
    parser.add_argument("-e", "--evalue", dest="evalue", required=True, type=float,
                        help="Maximum domain independent e-value to filter hmmscan output.", metavar="FLOAT")
    parser.add_argument("-v", "--validation", dest="validation", required=False, action='append',
                        help="Path to validation sample file. Parameter can be used repeatedly to include multiple samples.", metavar="FILE")
    parser.add_argument("-f", "--file", dest="file", nargs=2, action='append', default=[],
                        help="Config file variables to replace (e.g. --file PFAM2VEC path/to/pfam2vec.bin).", metavar="FILE")
    parser.add_argument("--log-dir", dest="log_dir", required=False, default='/tmp/bgc-pipeline-models',
                        help="Path to progress log directory (e.g. Tensorboard).", metavar="FILE")
    parser.add_argument("--log-file", dest="log_file", required=False,
                        help="Path to specific progress log file (e.g. Tensorboard).", metavar="FILE")
    parser.add_argument("--verbose", dest="verbose", required=False, default=2, type=int,
                        help="Verbosity level (0=none, 1=progress bar, 2=once per epoch).", metavar="INT")
    parser.add_argument(dest='samples', nargs='*',
                        help="Paths to training samples.", metavar="SAMPLES")
    options = parser.parse_args()

    with open(options.config, 'r') as fp:
        config = json.load(fp)

    progress_log_path = options.log_file or create_progress_log_path(options.log_dir, options.output)

    run_training(
        config=config,
        output_path=options.output,
        sample_paths=options.samples,
        validation_sample_paths=options.validation,
        evalue=options.evalue,
        progress_log_path=progress_log_path,
        files=options.file,
        verbose=options.verbose
    )