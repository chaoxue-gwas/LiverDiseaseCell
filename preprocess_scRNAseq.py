# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/04/07 21:18
# @Function: preprocess scRNA-seq dataset of liver.
# Liver Cell Atlas: https://www.livercellatlas.org/
# https://www.livercellatlas.org/download.php ; Liver Cell Atlas: Human: 1. All liver cells and 2. Visium spatial cells.
import gzip
import os
import argparse
import scanpy as sc
import scipy.io
import pandas as pd
from scipy import sparse

from para import RAW_DATA_DIR, PROJ_DATA_DIR
from util import make_dir


def read_10x_to_adata(matrix_dir,annot_h5ad_paths:[],cell_index_name='cell'):
    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    features_path = None
    for fname in ["features.tsv.gz", "genes.tsv.gz"]:
        try_path = os.path.join(matrix_dir, fname)
        if os.path.exists(try_path):
            features_path = try_path
            break
    if features_path is None:
        raise FileNotFoundError("Neither features.tsv.gz nor genes.tsv.gz found.")
    print("Reading matrix.mtx.gz...")
    with gzip.open(matrix_path, 'rt') as f:
        matrix = scipy.io.mmread(f).tocsc()
    print("Reading barcodes.tsv.gz...")
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].values
    print(f"Reading {os.path.basename(features_path)}...")
    with gzip.open(features_path, 'rt') as f:
        features = pd.read_csv(f, header=None, sep="\t")
    if features.shape[1] >= 2:
        gene_ids = features[0].values
        gene_names = features[1].values
    else:
        gene_ids = gene_names = features[0].values
    print("Constructing AnnData object...")
    adata = sc.AnnData(X=matrix.transpose())
    adata.var_names = gene_names
    adata.var["gene_ids"] = gene_ids
    adata.obs_names = barcodes
    for cell_annot_path,h5ad_out_path in annot_h5ad_paths:
        make_dir(os.path.dirname(h5ad_out_path))
        annotations = pd.read_csv(cell_annot_path)
        annotations.index=annotations[cell_index_name]
        common_cells = adata.obs_names.intersection(annotations.index)
        new_adata = adata[common_cells].copy()
        new_adata.obs = new_adata.obs.join(annotations, how='left')
        new_adata = new_adata[new_adata.obs['cluster'].notna()]
        new_adata.write(h5ad_out_path)
    # return adata


def main():
    # liver cells; 'humanAll'
    all_liver_cell_dir=f'{RAW_DATA_DIR}/scRNA-seq/rawData_human/countTable_human'
    annot_h5ads=[]
    for k in ['humanCD45neg','humanLymphoid','humanMyeloid']:
        annot_h5ads.append([f'{RAW_DATA_DIR}/scRNA-seq/annot_{k}.csv',f'{PROJ_DATA_DIR}/scRNA-seq/LCA_{k}.h5ad'])
    read_10x_to_adata(all_liver_cell_dir,annot_h5ads,'cell')
    # hepa cells
    hepa_cell_dir=f'{RAW_DATA_DIR}/scRNA-seq/rawData_humanVisium/countTable_humanVisium'
    hepa_cell_annot=f'{RAW_DATA_DIR}/scRNA-seq/annot_humanVisium.csv'
    hepa_cell_h5ad=f'{PROJ_DATA_DIR}/scRNA-seq/LCA_humanVisium.h5ad'
    # read_10x_to_adata(hepa_cell_dir,[hepa_cell_annot,hepa_cell_h5ad],'spot')


if __name__ == "__main__":
    main()


def main():
    pass


if __name__ == '__main__':
    main()
