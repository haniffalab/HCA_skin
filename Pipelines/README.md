# Pipelines
This Folder contains generalisable scripts which were used to analyse 10X(adult healthy, PS, AD and fetal skin) data, IHC scoring data and SS2 data.

## 01_Data_preprocessing
This directoery contains scripts used to import, build, qc and preprocess 10X cell ranger output data. It also includes ther harmony batch correction script and parameters used to batch correct our data. All datasets were processed in the same way.

## 02_leiden-cluster_consensus_logistic_regression_label_stability_prediction
This directory contains a custom logistic regression script to build label trasnfer models which may be used to test label stability at various continous dimension reduction parameters or upon selective removal of inputs within a categorical variable.T his script was used to compare analogous annotations between datasets to ensure accuracy and consistency of annotation and produces a heatmap weighted by probability of alignment derived from the binary assignment of all cells to a category by logistic regression. Written by IG. The script takes as input:

the analogous categorical variable in each dataset to be compared (e.g cell.labels)

## 03_AutoGeneS_implementation
A general implementation of the AutogeneS package

## 04_celltype_proportion_quasi_binom_neg_binom_statistical_testing
A general script which models the quasi binomial or negative binomial distribution of any categorical data (e.g celltype proportion).

## 05_VDJ_analysis_SCIRPY
A general script for VDJ analysis using the SCIRPY package.

## 06_cross_dataset_comparisons
A set of simple scripts to compute conserved markers between categoricals within a dataset

## 07_IHC_scoring statitistical models
A Simple script which takes a dataframe input containing observer quantified count data dervied from IHC and automatically selects paramateric or non-paramterics tests based on normality testing.

## 08_web_portal
Script used to generate the webportals available at (https://developmentcellatlas.ncl.ac.uk/)

## 09_interactive_heatmap
Script used to generate the interactive_heatmaps available at (https://developmentcellatlas.ncl.ac.uk/)

## 10_pseudotime_webportal
Script used to generate the pseudotime webportals available at (https://developmentcellatlas.ncl.ac.uk/)

## 11_add_dr_harmony_degs_annot
A general script which adds harmony corrected PC coordinates to a dataset

## 12_geneset_enrichment_annotation
This script was used to produce the chord plots. Chord plots are network visualisations representing neighbourhood structures of relationships between enriched genesets derived from GO.BP and other gene-function association databases. The script proceeds to cluster and annotate these pathways (please note that the cytoscape app has to be installed for this script to fully function). Written by IG. The script takes as input a ranked list of genes (e.g DEGs ranked by log fold change) and finds enriched pathways, it requires it's user to define:

The ranked genes as input (it is recommended that this list exceeds 50 genes)
The databases to acquire pathway-gene association data from (default is GO.BP)

## 13_hypergeometric_prediction_of_tfs
This script was used to compute predicted TFs which regulate genesets.(Output not shown)
