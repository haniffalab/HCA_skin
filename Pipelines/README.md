# Pipelines
This Folder contains generalisable scripts which were used to analyse 10X(adult healthy, PS, AD and fetal skin) data, IHC scoring data and SS2 data.

## 1_Data_preprocessing
This directoery contains scripts used to import, build, qc and preprocess 10X cell ranger output data. It also includes ther harmony batch correction script and parameters used to batch correct our data. All datasets were processed in the same way.

## 2_leiden-cluster_consensus_logistic_regression_label_stability_prediction
This directory contains a custom logistic regression script to build label trasnfer models which may be used to test label stability at various continous dimension reduction parameters or upon selective removal of inputs within a categorical variable.

