### Note!

**This repository provides data and examples that were used for development of DeepBGC and its evaluation with ClusterFinder and antiSMASH.**

**See https://github.com/Merck/deepbgc for the DeepBGC tool.**

### Note!

# DeepBGC development & evaluation code

## Reproducing data

Reproduction and storage of data files is managed using [DVC](https://github.com/iterative/dvc) (development version `0.22.0`). 
Each data file has a `.dvc` history file that contains the command that was used to generate the output along with md5 hashes of its dependencies.

## Installation

- Install python 3, ideally using conda
- Run `pip install -r requirements.txt` to download DVC and other requirements

## Downloading a file

- Run the AWS config script to generate temporary AWS credentials in ~/.aws/credentials:
  - `generate-aws-config --account lab --insecure`
- Run `dvc pull data/path/to/file.dvc` to download required file.

## High-level overview

### Main folders

- [bgc_detection/](bgc_detection/) all the code
- [data/](notebooks/) all the data
    - [bacteria/](data/bacteria) 3k reference bacteria
        - [candidates/](data/bacteria/candidates) novel detected BGC candidates
    - [clusterfinder/](data/clusterfinder) ClusterFinder (Cimermancic et al.) datasets
    - [evaluation/](data/evaluation) Cross-validation, Leave-Class-Out and Bootstrap evaluation
    - [features/](data/features) Pfam2vec and other protein domain features
    - [figures/](data/figures) Paper figures
    - [mibig/](data/mibig) MIBiG BGC database samples
    - [models/](data/models) Model configurations and trained models
    - [pfam/](data/pfam) Pfam repository files
    - [training/](data/training) Negative and positive training data and t-SNE visualizations
- [notebooks/](notebooks/) Jupyter (iPython) notebooks


### Training a model

- Define a JSON config file, see [data/models/config](data/models/config) for reference.
- Run [bgc_detection/run_training.py](bgc_detection/run_training.py) with given config and path to training data. See DVC files in [data/models/trained](data/models/trained) for reference.
- Trained model will be presented as Python pickle file. 

### Predicting using trained model

- Prepare a protein FASTA file, e.g. using [Prodigal](https://github.com/hyattpd/Prodigal) 
(see [data/bacteria/proteins.dvc](data/bacteria/proteins.dvc) for reference) or extract it from an annotated GenBank file using [bgc_detection/preprocessing/proteins2fasta.py](bgc_detection/preprocessing/proteins2fasta.py).
- Detect protein domains using Hmmscan (see [data/bacteria/domtbl.dvc](data/bacteria/domtbl.dvc) for reference)
- Convert the Hmmscan domtbl file into a Domain CSV file using [bgc_detection/preprocessing/domtbl2csv.py](bgc_detection/preprocessing/domtbl2csv.py) 
(see [data/bacteria/domains.dvc](data/bacteria/domains.dvc) for reference)
- Predict BGC domain-level probability using [bgc_detection/run_prediction.py](bgc_detection/run_prediction.py) 
(see [data/bacteria/prediction/128lstm-100pfamdim-8pfamiter-posweighted-neg-10k.dvc](data/bacteria/prediction/128lstm-100pfamdim-8pfamiter-posweighted-neg-10k.dvc) for reference)
- Threshold and merge domain-level predictions into a BGC candidate CSV file using 
[bgc_detection/candidates/threshold_candidates.py](bgc_detection/candidates/threshold_candidates.py) 
(see [data/bacteria/candidates/128lstm-100pfamdim-8pfamiter-posweighted-neg-10k-fpr2/candidates.csv.dvc] for reference)

### Bootstrap validation on 9 Fully-annotated genomes

See [notebooks/LabelledContigBootstrap.ipynb](notebooks/LabelledContigBootstrap.ipynb).

### Leave Class Out validation and Cross validation

See [data/evaluation/lco-neg-10k](data/evaluation/lco-neg-10k) (TODO).

See [data/evaluation/cv-10fold-neg-10k](data/evaluation/cv-10fold-neg-10k) (TODO).

### Random Forest classification

See [notebooks/CandidateClassification.ipynb](notebooks/CandidateClassification.ipynb) and [notebooks/CandidateActivityClassification.ipynb](notebooks/CandidateActivityClassification.ipynb)

### Novel BGC candidates generation

See [notebooks/NovelCandidates.ipynb](notebooks/NovelCandidates.ipynb).
