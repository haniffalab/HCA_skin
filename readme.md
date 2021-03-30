# Developmental cell programs are co-opted in inflammatory skin diseases

This repository contains code used for the analysis of the Human Cell Atlas Skin paper published in Science Jan 2021 DOI: 10.1126/science.aba6500

## Data
Datasets analysed include: 10X data (Healthy adult Skin, PS skin, AD skin, fetal Skin), IHC for healthy and diseased skin, SS2 data, and healthy skin BCR-/TCR-enriched VDJ data.

To quickly browse our data, our online webportal resource may be accessed via (https://developmentcellatlas.ncl.ac.uk/)

Code is also available on Zenodo 10.5281/zenodo.4249674

Data is available as a h5ad file at 10.5281/zenodo.4569496

Raw data is available on ArrayExpress www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8142

## Analysis workflows
Cellranger count matrix files for Healthy, AD and PS skin 10X data were loaded into one object as described in https://github.com/haniffalab/FCA_liver for downstream analysis. 

### Pipelines
Generalisable scripts are saved in the 'pipelines' directory. Please refer to pipelines/readme file for further information on the methods used for each pipeline. The figure panels created through use of each pipeline are detailed in readme file.

### Celltype compartment analysis
Custom scripts for specific celltype compartment analysis that do not fall under the "generalisable" category are available here. 

All scripts written by IG and VEGP and GR unless otherwise stated.

See [Human Cell Atlas](https://www.humancellatlas.org) for more details.
