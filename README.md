<div align="center">    
 
# Investigating the ability of deep learning-based structure prediction to extrapolate and/or enrich the set of antibody CDR canonical forms

</div>

This repository contains code used to cluster CDR loops of predicted structures from paired sequences in Observed Antibody Space (OAS). In addition to the code we provide the sequence motifs of density based clusters found in our analyses.

## Sequence motifs of density based clusters
Sequence motifs and density based clusters are given as SVG files. Folders are structured by CDR and length.

## Clustering pipeline

The git folder structure must be maintained to ensure outputs are placed into correct folders for subsequent steps of the pipeline to work.

### Download the necessary structures and resource files.

- Predicted structures can be downloaded from: https://doi.org/10.5281/zenodo.10280180
- The .tar.gz files containing the predicted structures should be downloaded to ``` input/oas_structures ```. A small sample of around 100K should be fine to run the code, but this may produce different results (clusters and logo plots) in comparison to when running on all 1.5M.
- SAbDab structures and resource directories can be downloaded here: https://doi.org/10.5281/zenodo.10554753
- Both SAbDab reference structures (```imgt_sabdab```) and the split structures (```sabdab_chains_renumb```) which have been curated into single chains and saved as individual pdb files per chain, are needed to run the code. These directories should be placed in the input folder (see placeholders). 
- The resources (```resources```) folder contains the reference csv files, and RDS sequence files for downstream analysis. This should be at the top level (see placeholder).




The files can be unpacked following the download by running the following command in the oas_structures folder after all files have been downloaded:

```

# Unpack the tar balls.
for file in *.tar.gz; do
     echo "Extracting $file..."
     tar -xzf "$file" -C .
done
echo "Extraction complete."

# Unpack all sub folders to top level.
find "." -mindepth 2 -type f -name "*.pdb" -exec mv {} "." 
```

NB: If you want to run a custom analysis with a specific loop - ensure that only that loop is uncommented in the ```_01_align_pipeline.py``` file.

Files are run in the following order.
1. Modify ```_01_align_pipeline.py``` to select the loop of interest. E.g uncomment ```('L', 3, 10)```. This script aligns all relevant loops to a template. These will be output to output/aligned/ directory in sub directories names by the relevant cdr and length.
2. ```_02_rmsd_pairwise.py``` calculates the pairwise RMSD between all aligned loops and saves the results in a single long csv file.
3. ```_03_run_mds.R``` converts this to a pairwise matrix and runs multi-dimensional scaling on the matrix to produce the MDS corordinates.
4. ```_04_dbscan.py``` carries out density based clustering on the pairwise matrix and outputs plots for different values of K as well as the cluster assignment for each value of K.
5. ```_05_mds_vis.R``` visualises the MDS coordinates, creates logo plots and overlays the MDS plot with PyIgClassify2 assignments.

These can be run with the script:
run_pipeline.sh
