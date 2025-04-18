# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/04/07 11:29
# @Function: run scDRS

import os
import pickle
import re

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scdrs
import scanpy as sc
import scipy.sparse as sparse
import urllib3
from anndata import read_h5ad
from scdrs.util import plot_heatmap, small_squares
from statsmodels.stats.multitest import multipletests

from para import PROJ_DATA_DIR, PROJ_RESULT_DIR
from scdrs_util import plot_group_stats
from util import make_dir, log, go_enrich_plot


urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
os.environ['HTTP_PROXY'] = 'socks5://127.0.0.1:7898'
os.environ['HTTPS_PROXY'] = 'socks5://127.0.0.1:7898'


def load_h5ad(
    h5ad_file: str, flag_filter_data: bool = False, flag_raw_count: bool = True
) -> anndata.AnnData:
    """Load h5ad file and optionally filter out cells and perform normalization.

    Parameters
    ----------
    h5ad_file : str
        Path to h5ad file
    flag_filter_data : bool
        If True, filter out cells with

        - sc.pp.filter_cells(adata, min_genes=250)
        - sc.pp.filter_genes(adata, min_cells=50)
    flag_raw_count : bool
        If True, perform size-factor normalization and log1p transformation.

    Returns
    -------
    adata : anndata.AnnData
        Single-cell data.
    """
    adata = read_h5ad(h5ad_file)
    # check inputs (1) no NaN in adata.X (2) all adata.X >= 0
    if np.isnan(adata.X.sum()):
        raise ValueError(
            "h5ad expression matrix should not contain NaN. Please impute them beforehand."
        )
    if (adata.X < 0).sum() > 0:
        raise ValueError(
            "h5ad expression matrix should not contain negative values. "
            "This is because in the preprocessing step, "
            "scDRS models the gene-level log mean-variance relationship. "
            "See scdrs.pp.compute_stats for details."
        )

    if flag_filter_data:
        sc.pp.filter_cells(adata, min_genes=250)
        sc.pp.filter_genes(adata, min_cells=50)
    if flag_raw_count:
        if sparse.issparse(adata.X):
            adata.X = adata.X.astype(np.float64)
        else:
            adata.X = np.array(adata.X, dtype=np.float64)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
    return adata


def rename_var_conflicts(adata, obs_keys, suffix='_gene'):
    var_names = adata.var_names.tolist()
    rename_map = {}

    for key in obs_keys:
        if key in var_names:
            new_name = f"{key}{suffix}"
            while new_name in var_names or new_name in rename_map.values():
                new_name += '_'
            rename_map[key] = new_name

    if rename_map:
        for old, new in rename_map.items():
            print(f" - {old} → {new}")
        adata.var_names = [rename_map.get(name, name) for name in adata.var_names]
    return rename_map


def run(sc_data_key):
    scdrs_result_pyd=f'{PROJ_RESULT_DIR}/scdrs/assoc_cell/{sc_data_key}/{sc_data_key}.pyd'
    H5AD_FILE=f'{PROJ_DATA_DIR}/scRNA-seq/{sc_data_key}.h5ad'
    GS_FILE=f'{PROJ_RESULT_DIR}/scdrs/gs/processed_geneset.gs'
    adata = load_h5ad(H5AD_FILE, flag_filter_data=True, flag_raw_count=True)
    df_gs = scdrs.util.load_gs(GS_FILE,to_intersect=adata.var_names)
    scdrs.preprocess(adata)
    res_map={}
    if os.path.isfile(scdrs_result_pyd):
        with open(scdrs_result_pyd, 'rb') as f:
            res_map:dict=pickle.load(f)
            log(f'load exsited {len(res_map)} phenos')
    exsited_phenos=sorted(res_map.keys())
    for pheno in sorted(df_gs.keys()):
        # if pheno in exsited_phenos:
        #     continue
        if pheno not in ['PBC']:
            continue
        log(f'start run: {pheno}')
        gene_list = df_gs[pheno][0]
        gene_weight = df_gs[pheno][1]
        df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight,
                                  ctrl_match_key="mean_var",
                                  n_ctrl=1000,
                                  weight_opt="vs",
                                  return_ctrl_raw_score=False,
                                  return_ctrl_norm_score=True,
                                  verbose=False,
                                  )
        res_map[pheno]=df_res
    make_dir(os.path.dirname(scdrs_result_pyd))
    with open(scdrs_result_pyd, 'wb') as f:
        pickle.dump(res_map, f)
        log(f'save to: {scdrs_result_pyd}')

def analyze(sc_data_key,cell_annot_col):
    res_dir=f'{PROJ_RESULT_DIR}/scdrs/assoc_cell/{sc_data_key}'
    scdrs_result_pyd=f'{res_dir}/{sc_data_key}.pyd'
    scdrs_group_pyd=f'{res_dir}/{sc_data_key}.cell_type.pyd'
    H5AD_FILE=f'{PROJ_DATA_DIR}/scRNA-seq/{sc_data_key}.h5ad'
    adata = load_h5ad(H5AD_FILE, flag_filter_data=True, flag_raw_count=True)
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata,n_neighbors = 15, n_pcs = 20)
    scdrs.preprocess(adata)
    log(f'load sc data')
    res_ct_map={}
    if os.path.isfile(scdrs_group_pyd):
        with open(scdrs_group_pyd, 'rb') as f:
            res_ct_map:dict=pickle.load(f)
            log(f'load exsited {len(res_ct_map)} phenos')
    exsited_phenos=sorted(res_ct_map.keys())
    with open(scdrs_result_pyd, 'rb') as f:
        res_map:dict=pickle.load(f)
        log(f'load scdrs results')
        for pheno in sorted(res_map.keys()):
            # if pheno in exsited_phenos:
            #     continue
            if pheno not in ['PBC']:
                continue
            log(f'start analyze {pheno}')
            # code for analysis
            df_stats = scdrs.method.downstream_group_analysis(
                adata=adata,
                df_full_score=res_map[pheno],
                group_cols=[cell_annot_col],
            )[cell_annot_col]
            res_ct_map[pheno]=df_stats
            # print(df_stats)
    make_dir(os.path.dirname(scdrs_group_pyd))
    with open(scdrs_group_pyd, 'wb') as f:
        pickle.dump(res_ct_map, f)
        log(f'save to: {scdrs_group_pyd}')

def plot_legend(cell_types,palette):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               label=ct, markerfacecolor=palette[i], markersize=10)
        for i, ct in enumerate(cell_types)
    ]
    fig, ax = plt.subplots(figsize=(4, len(cell_types) * 0.3 + 1))
    ax.axis('off')
    ax.legend(handles=legend_elements, loc='center', frameon=False, ncol=1)
    plt.tight_layout()
    plt.show()

def write_result_table(sc_data_key):
    res_dir=f'{PROJ_RESULT_DIR}/scdrs/assoc_cell/{sc_data_key}'
    out_excel_path=f'{PROJ_RESULT_DIR}/scdrs/analysis/assoc_cell_type/{sc_data_key}.xlsx'
    scdrs_group_pyd=f'{res_dir}/{sc_data_key}.cell_type.pyd'
    if not os.path.isfile(scdrs_group_pyd):
        raise Exception('no result file')
    with open(scdrs_group_pyd, 'rb') as f:
        res_ct_map:dict=pickle.load(f)
        make_dir(os.path.dirname(out_excel_path))
        with pd.ExcelWriter(out_excel_path) as writer:
            for phenotype, df in res_ct_map.items():
                sheet_name = phenotype
                df.index.name='cell_type'
                df.to_excel(writer, sheet_name=sheet_name)
                log(f'save to {out_excel_path}')



def visualize_common(sc_data_key,cell_annot_col):
    plot_phenos=['HBV','HCV','PBC','PSC','NAFLD','SJ','TC','TG'] #,'TC','TG'
    res_dir=f'{PROJ_RESULT_DIR}/scdrs/assoc_cell/{sc_data_key}'
    scdrs_result_pyd=f'{res_dir}/{sc_data_key}.pyd'
    scdrs_group_pyd=f'{res_dir}/{sc_data_key}.cell_type.pyd'
    H5AD_FILE=f'{PROJ_DATA_DIR}/scRNA-seq/{sc_data_key}.h5ad'
    adata = load_h5ad(H5AD_FILE, flag_filter_data=True, flag_raw_count=True)
    adata.obsm['X_umap'] = adata.obs[['UMAP_1', 'UMAP_2']].to_numpy()
    with open(scdrs_result_pyd, 'rb') as f:
        dict_score:dict=pickle.load(f)
        for trait in dict_score:
            adata.obs[trait] = dict_score[trait]["norm_score"]
        sc.set_figure_params(figsize=[6, 6], dpi=150)
        sc.pl.umap(
            adata,
            color=cell_annot_col,
            s=5,
            legend_loc = None
        )
        cell_types = adata.obs[cell_annot_col].astype('category').cat.categories
        cell_count=adata.obs[cell_annot_col].value_counts().to_dict()
        colors = adata.uns[f'{cell_annot_col}_colors']
        cts=[f'{c} ({cell_count[c]})' for c in cell_types]
        plot_legend(cts,colors)
        rename_var_conflicts(adata,sorted(dict_score.keys()))
        sc.set_figure_params(figsize=[3, 3], dpi=150)
        sc.pl.umap(
            adata,
            color=plot_phenos,
            color_map="RdBu_r",
            vmin=0,
            vmax=6,
            s=20,
            ncols=4,
        )
        plt.show()
    with open(scdrs_group_pyd, 'rb') as f:
        res_ct_map:dict=pickle.load(f)
        clean_map={p:res_ct_map[p] for p in plot_phenos}
        # print(res_ct_map)
        fig,ax=plot_group_stats(
            dict_df_stats=clean_map,
            plot_kws={
                "cb_vmax": 1,
                "cb_fraction": 0.12,
                "cb_n_bin": 10,
            }
        )
        plt.grid(False)
        plt.show()
    pass

class HeteroAnalysis:
    def __init__(self):
        self.plot_phenos = ['HBV', 'HCV', 'PBC', 'PSC', 'NAFLD', 'SJ', 'TC', 'TG']  # ,'TC','TG'
        self.out_common_dir = f'{PROJ_RESULT_DIR}/scdrs/analysis/hetero_cell_common'
        self.out_function_dir= f'{PROJ_RESULT_DIR}/scdrs/analysis/hetero_cell_function'

    def ab_function(self,adata_ca1,eff_groups,out_prefix,top_n=200):
        sc.set_figure_params(figsize=[6, 6], dpi=150)
        adata_ca1 = adata_ca1.copy() #[adata_ca1.obs['merged_groups'].isin([eff_group,ref_group])]
        sc.tl.rank_genes_groups(adata_ca1, groupby='merged_groups', groups=eff_groups,method='wilcoxon')
        result = adata_ca1.uns['rank_genes_groups']
        make_dir(os.path.dirname(out_prefix))
        for clu in eff_groups:
            degs = pd.DataFrame({
                'gene': result['names'][clu],
                'pval_adj': result['pvals_adj'][clu],
                'score': result['scores'][clu],
            })
            df_sorted = degs.sort_values(by='score', ascending=False)
            ga_genes = df_sorted.head(top_n)['gene'].values
            ge=[ga_genes,f'{out_prefix}.{clu}.enrich.jpg',clu]
            go_enrich_plot([str(g) for g in ge[0]],ge[1],title=ge[2],anno_dbs=['KEGG','GO:BP'],plot_max_term=10)
            df_sorted.head(top_n).to_excel(f'{out_prefix}.{clu}.dif_genes.top{top_n}.xlsx',index=False)


    def analysis_common(self,sc_data_key,plot_phenos,merge_cluster_arr:dict,select_celltype,out_prefix,plot_size=10):
        res_dir = f'{PROJ_RESULT_DIR}/scdrs/assoc_cell/{sc_data_key}'
        scdrs_result_pyd = f'{res_dir}/{sc_data_key}.pyd'
        H5AD_FILE = f'{PROJ_DATA_DIR}/scRNA-seq/{sc_data_key}.h5ad'
        with open(scdrs_result_pyd, 'rb') as f:
            dict_score: dict = pickle.load(f)
            adata = load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
            rename_var_conflicts(adata, sorted(dict_score.keys()))
            adata_ca1 = adata[adata.obs[cell_annot_col].isin([select_celltype])].copy()
            # re-cluster
            sc.pp.filter_cells(adata_ca1, min_genes=0)
            sc.pp.filter_genes(adata_ca1, min_cells=1)
            sc.pp.normalize_total(adata_ca1, target_sum=1e4)
            sc.pp.log1p(adata_ca1)
            sc.pp.highly_variable_genes(adata_ca1, n_top_genes=2000)
            adata_ca1 = adata_ca1[:, adata_ca1.var.highly_variable]
            sc.pp.scale(adata_ca1, max_value=10)
            sc.tl.pca(adata_ca1, svd_solver="arpack")
            sc.pp.neighbors(adata_ca1, n_neighbors=10, n_pcs=40)
            sc.tl.leiden(adata_ca1, resolution=0.5)
            sc.tl.umap(adata_ca1, n_components=2)
            # assign scDRS score
            for trait in dict_score:
                adata_ca1.obs[trait] = dict_score[trait]["norm_score"]
            sc.set_figure_params(figsize=[3, 3], dpi=150)
            sc.pl.umap(
                adata_ca1,
                color=plot_phenos,
                color_map="RdBu_r",
                vmin=0,
                vmax=6,
                s=20,
                ncols=4,
            )
            if merge_cluster_arr is None:
                adata_ca1.obs['merged_groups']=adata_ca1.obs['leiden'].astype('category')
                annot_clusters=adata_ca1.obs['merged_groups'].unique().tolist()
            else:
                merge_cluster = {}
                for k in merge_cluster_arr.keys():
                    for v in merge_cluster_arr[k]:
                        merge_cluster[v] = k
                for clu in adata_ca1.obs['leiden'].unique():
                    if clu not in merge_cluster.keys():
                        merge_cluster[clu]=''
                adata_ca1.obs['merged_groups'] = adata_ca1.obs['leiden'].map(merge_cluster).astype('category')
                annot_clusters=sorted(merge_cluster_arr.keys())
            sc.set_figure_params(figsize=[3, 3], dpi=150)
            sc.pl.umap(adata_ca1, color='leiden', s=plot_size, legend_loc='on data')
            sc.set_figure_params(figsize=[3, 3], dpi=150)
            sc.pl.umap(adata_ca1, color='merged_groups', s=plot_size, legend_loc='on data')
            # function analysis
            self.ab_function(adata_ca1,annot_clusters,out_prefix)

    def plot_hereto_cell(self,sc_data_key,select_celltypes,out_prefix):
        res_dir = f'{PROJ_RESULT_DIR}/scdrs/assoc_cell/{sc_data_key}'
        scdrs_result_pyd = f'{res_dir}/{sc_data_key}.pyd'
        H5AD_FILE = f'{PROJ_DATA_DIR}/scRNA-seq/{sc_data_key}.h5ad'
        make_dir(os.path.dirname(out_prefix))
        for select_celltype in select_celltypes:
            with open(scdrs_result_pyd, 'rb') as f:
                dict_score: dict = pickle.load(f)
                adata = load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
                rename_var_conflicts(adata, sorted(dict_score.keys()))
                adata_ca1 = adata[adata.obs[cell_annot_col].isin([select_celltype])].copy()
                # re-cluster
                sc.pp.filter_cells(adata_ca1, min_genes=0)
                sc.pp.filter_genes(adata_ca1, min_cells=1)
                sc.pp.normalize_total(adata_ca1, target_sum=1e4)
                sc.pp.log1p(adata_ca1)
                sc.pp.highly_variable_genes(adata_ca1, n_top_genes=2000)
                adata_ca1 = adata_ca1[:, adata_ca1.var.highly_variable]
                sc.pp.scale(adata_ca1, max_value=10)
                sc.tl.pca(adata_ca1, svd_solver="arpack")
                sc.pp.neighbors(adata_ca1, n_neighbors=10, n_pcs=40)
                sc.tl.leiden(adata_ca1, resolution=0.5)
                sc.tl.umap(adata_ca1, n_components=2)
                # assign scDRS score
                for trait in dict_score:
                    adata_ca1.obs[trait] = dict_score[trait]["norm_score"]
                sc.set_figure_params(figsize=[3, 3], dpi=150)
                ct_str=re.sub(r'[^a-zA-Z0-9]', '_', select_celltype)
                sc.settings.figdir = os.path.dirname(out_prefix)
                sc.pl.umap(
                    adata_ca1,
                    color=self.plot_phenos+['leiden'],
                    color_map="RdBu_r",
                    vmin=0,
                    vmax=6,
                    s=20,
                    ncols=4,
                    show=False,
                    legend_loc='on data',
                )
                plt.suptitle(f'Disease score of {select_celltype}', fontsize=17)
                plt.tight_layout(rect=[0, 0, 1, 1])
                plt.savefig(f'{out_prefix}.{ct_str}.png', dpi=300)
                plt.close()


    def analysis_all_hetero_celltypes(self):
        hereto_params={
            'LCA_humanCD45neg':['Central Vein Endothelial cells','Cholangiocytes','Fibroblasts','Hepatocytes',
                                'LSECs','Portal Vein Endothelial cells'],
            'LCA_humanLymphoid':['B cells','Circulating TEM','Cytotoxic CD8+','Gd T cells','RM CD8+ T cells','pDCs'],
            'LCA_humanMyeloid':['MoMac1','Pre-moKCs and moKCs','cDC1s','cDC2s','resKCs']
        }
        out_dir=f'{self.out_common_dir}'
        for sc_data_key in hereto_params.keys():
            self.plot_hereto_cell(sc_data_key,hereto_params[sc_data_key],f'{out_dir}/{sc_data_key}')


    def analysis_1(self):
        merge_cluster={'GroupA':['0','1'],'GroupB':['2','3','4','5']}
        select_ct="LSECs"
        sc_data_key='LCA_humanCD45neg'
        out_prefix=f'{self.out_function_dir}/{sc_data_key}.{select_ct}/{select_ct}'
        self.analysis_common(sc_data_key,self.plot_phenos,merge_cluster,select_ct,out_prefix)

    def analysis_2(self):
        merge_cluster={'GroupA':['0','1','2','5'],'GroupB':['3','4','6','7']}
        select_ct="Cholangiocytes"
        sc_data_key='LCA_humanCD45neg'
        plot_phenos=['PBC', 'PSC', 'TC']
        out_prefix=f'{self.out_function_dir}/{sc_data_key}.{select_ct}/{select_ct}'
        self.analysis_common(sc_data_key,plot_phenos,merge_cluster,select_ct,out_prefix,20)

    def analysis_3(self):
        merge_cluster={'1':['1'],'2':['2']}
        # merge_cluster=None
        select_ct="RM CD8+ T cells"
        sc_data_key='LCA_humanLymphoid'
        plot_phenos=['HBV','HCV','PBC', 'PSC']
        out_prefix=f'{self.out_function_dir}/{sc_data_key}.{select_ct}/{select_ct}'
        self.analysis_common(sc_data_key,plot_phenos,merge_cluster,select_ct,out_prefix,20)



if __name__ == '__main__':
    # sc_data_key = 'LCA_human_all_liver_cell'  # LCA_human_spatial_cell
    cell_annot_col = 'annot'  # annot
    sc_keys=['humanCD45neg','humanLymphoid','humanMyeloid']
    for k in sc_keys:
        log(f'start run {k}')
        sc_data_key=f'LCA_{k}'
        # run(sc_data_key)
        # analyze(sc_data_key,cell_annot_col)
        # write_result_table(sc_data_key)
        # visualize_common(sc_data_key, cell_annot_col)
    # hero
    ha=HeteroAnalysis()
    # ha.analysis_all_hetero_celltypes()
    ha.analysis_1()
    ha.analysis_2()
    ha.analysis_3()