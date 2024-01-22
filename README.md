<div align="center">    
 
# Investigating the ability of deep learning-based structure prediction to extrapolate and/or enrich the set of antibody CDR canonical forms

</div>

This repository contains code used to cluster CDR loops of predicted structures from paired sequences in Observed Antibody Space (OAS). In addition to the code we provide the sequence motifs of density based clusters found in our analyses.

## Sequence motifs of density based clusters
Sequence motifs and density based clusters are given as SVG files. Folders are structured by CDR and length.

## Clustering pipeline
- Predicted structures can be downloaded from: https://zenodo.org/records/10280181
- To facilitate analysis of SAbDab structures, they have been curated into single chains and saved as individual pdb files per chain. These can be downloaded here:
- For assignment of PyIgClassify2 canonical cluster information the file: pyig_cdr_data.txt.gz can be downloaded from: http://dunbrack2.fccc.edu/PyIgClassify2/ this must be unzipped and placed in resources/

The above folder structure must be maintained to ensure outputs are placed into correct folders for subsequent steps of the pipeline to work.
To run the clustering pipeline, ensure ABB2 structures (a small sample is fine, but may produce different results to when running on all) are downloaded and unpacked into: input/oas_structures/. Ensure the curated SAbDab structures have been saved to: input/sabdab_chains_renumb/

All output files will be saved with the folder structure in a location that allows the next step to pick up and run all files in the relevant folder.

Files are run in the following order.
1. Modify _01_align_pipeline.py to select the loop of interest. E.g uncomment ('L', 3, 10). Run this file to produce the aligned loops. These will be output to output/aligned/ followed by the relevant cdr and length.
2. Run the _02_rmsd_pairwise.py file. This calculates the pairwise RMSD between all aligned loops and saves the results in a single long csv file.
3. _03_run_mds.R converts this a pairwise matrix and runs multi-dimensional scaling on the matrix to produce the MDS corordinates.
4. _04_dbscan.py carries out density based clustering on the pairwise matrix and outputs plots for different values of K as well as the cluster assignment for each value of K.
5. _05_mds_vis.R visualised the MDS coordinates, creates logo plots and overlays the MDS plot with PyIgClassify2 assignments.

These can be run with the script:
run_pipeline.sh
