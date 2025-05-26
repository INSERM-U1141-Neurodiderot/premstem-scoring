# Systematic outcome scoring for enhancing in vivo neurotherapeutic testing applied to perinatal brain injury

Official code for the analyses performed in the paper "Systematic outcome scoring for enhancing in vivo neurotherapeutic testing applied to perinatal brain injury". This code is under a license approved by the open source initiative.

[![Python Version](https://img.shields.io/badge/python-3.7-pink)](https://img.shields.io/badge/python-3.8%7C3.9-pink) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 1. System requirements

**OS:** Linux Debian
**Dependencies**: file *requirements.txt*

## 2. Citation

```
@Unpublished{bokobza2025systematic,
  title={Systematic outcome scoring for enhancing in vivo neurotherapeutic testing applied to perinatal brain injury},
  author={Bokobza, Cindy and Réda, Clémence and Nair, Syam et al.},
  note={Under review},
  year={2025},
}
```

## 3. Installation guide (estimated install time: ~1min)

It is strongly advised to install the [Conda](https://docs.anaconda.com/free/miniconda/miniconda-install/) tool to create a virtual environment, and [Pip](https://pip.pypa.io/en/stable/installation/) for installing dependencies:

```bash
conda create --name scoring python=3.6 -y
conda activate scoring
python3 -m pip install -r requirements.txt
```

## 4. Demo (estimated runtime: <2h)

### 4.a. Download the data

Download the data as described in the paper.

Filename             |    Administration protocol (time and mode)
---------------------|---------------------------------------------
P7 (inas P5)         | P7 / INAS / P5
results_user_run409  | P7 / IV  / P5
results_user_run412  | P12 / INAS / P10
results_user_run413  | P12 / IV  / P10
results_user_run415  | P22 / IV  / P20
results_user_run416  | P22 / INAS / P20

### 4.b. Compute Characteristic Direction (CD) signatures and the cosine scores

```bash
conda activate scoring
## Apply the pipeline described in the paper for each batch
python3 -m ranking_conditions "P7 (inas P5)"
python3 -m ranking_conditions "results_user_run409"
python3 -m ranking_conditions "results_user_run412"
python3 -m ranking_conditions "results_user_run413"
python3 -m ranking_conditions "results_user_run415"
python3 -m ranking_conditions "results_user_run416"
```

After launching the command
```bash
python3 -m ranking_conditions "<batch_name>"
```
a folder named "<batch_name>" will be created, containing (1) a file named "CD\_signatures.csv" with the signatures for each condition in the batch (PBS, Dose (low, medium, high)), (2) a file named "dose\_ranking.csv", (3) 4 individual files containing each signature (PBS, Dose1, Dose2, Dose3). The concatenation of each "dose\_ranking.csv" gives the ranking table with cosine scores.

### 4.c. Compare the ranking scores (cosine similarity) with the size of the gene support 

The support for a cosine score on two signatures is the set of genes which are present in both signatures.

```bash
conda activate scoring
python3 -m correlation_score_geneset
```

This command returns a figure "correlation_score_intersections.png" plotting cosine scores against the size of the gene support, and computes the R2 metric between those two measures (R2=0,13) and the Pearson's R (r=-0.36, p-value=0.14).

### 4.d. Compare CD signatures with DESeq-based signatures (histogram)

```bash
conda activate scoring
python3 -m compare_analyses "" CD ## plot figures for CD signatures
python3 -m compare_analyses "" DESeq2 ## plot figures for DESeq2 signatures
```
