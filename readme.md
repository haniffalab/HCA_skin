# HCA Skin single cell analysis

This repository contains code used for the analysis of the Human Cell Atlas Skin samples *(Fletcher, Vegh, Reynolds et al., 2019, manuscript submitted).*

All files in the 'R' directory were sourced prior to analysis. Unprocessed count files were loaded into one object as described in the 'data' directory, then cells were subset into four major groups and analysed. Finally, cell type annotations and raw counts were exported and uploaded to ArrayExpress. Please see `sessioninfo.txt` for the version of R packages used.

The web portal (https://developmentcellatlas.ncl.ac.uk/access/hca_skin_access) files were made using the single cell analysis bundle ([v1.0.0](https://github.com/haniffalab/scRNA-seq_analysis/releases/tag/v1.0.0), https://github.com/haniffalab/scRNA-seq_analysis), set up on the School of Computing Science HPC Facility, Newcastle University.

See [Human Cell Atlas](https://www.humancellatlas.org) for more details.
