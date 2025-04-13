##  Analysis Code for Liver Disease-Associated Cells
* ```preprocess_GWAS.py``` Prepares GWAS summary statistics for MAGMA, runs MAGMA, and extracts the top 1000 associated genes with their Z-scores for use in scDRS.
* ```preprocess_scRNAseq.py``` Converts human liver single-cell data from the LCA dataset into .h5ad format for use in scDRS.
* ```run_scDRS.py``` Runs scDRS for analysis and visualization.
* ```para.py``` Stores parameters and paths to resource files.
* ```util.py``` Utility functions.
* ```scdrs_util.py``` Utility functions from scDRS v1.0.3.
