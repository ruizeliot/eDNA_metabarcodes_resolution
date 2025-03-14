# eDNA_metabarcodes_resolution

Scripts, functions and graphs accompanying the whole R project with data and outputs deposited on: https://doi.org/10.6084/m9.figshare.26014840

All R functions are commented and provided with a manual (see "man" folder). 

These files can be viewed in Rstudio (click on "Preview" after opening the file) under the same format than the help of a function within a package. 

They are also provided with small and reproducible examples to show how the function works.

The program VSEARCH (https://github.com/torognes/vsearch?tab=readme-ov-file) needs to be installed (and its installation path indicated) to run the functions "vsearch_pairwise_similarity", "intra_taxa_similarity" and "multiple_resolution_analyses".

The programs MOTHUR (https://mothur.org/wiki/installation/) and SWARM (https://github.com/torognes/swarm/releases) need to be installed (and their installation path indicated) to run the function "cluster_metabarcodes_mothur_swarm".
