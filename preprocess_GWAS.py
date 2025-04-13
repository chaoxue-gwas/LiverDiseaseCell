# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/04/07 11:30
# @Function: preprocess GWAS summary data
import os.path

import pandas as pd

from para import RAW_DATA_DIR, PROJ_DATA_DIR, magma_gene_annot, PROJ_RESULT_DIR, magma_ref_genos, snp_rsid_db_hg19
from util import log, make_dir, LOCAL_DIR, batch_shell_plus


def make_magam_gwas_format():
    log(f'start')
    sc_dir=f'{RAW_DATA_DIR}/gwas'
    out_dir=f'{PROJ_DATA_DIR}/gwas'
    make_dir(out_dir)
    meta_path=f'{sc_dir}/gwas.meta.xlsx'
    mdf=pd.read_excel(meta_path)
    snp_map = pd.read_table(snp_rsid_db_hg19, sep='\t',dtype=str)
    snp_map.columns=['CHR','BP','rsid']
    log(f'load {snp_map.shape[0]} snps from common snp db')
    for i in mdf.index:
        abbr=mdf.loc[i,'abbr']
        out_path=f'{out_dir}/{abbr}.magma.gwas.tsv'
        if os.path.isfile(out_path):
            log(f'{abbr} existed, skip')
            continue
        log(f'start {abbr}')
        sep=str(mdf.loc[i,'sep']).strip('"')
        df=pd.read_table(f'{sc_dir}/{mdf.loc[i,"path"]}',dtype=str,sep=sep)
        cols=['snp','chr','bp','p','n_eff']
        final_cols=['SNP','CHR','BP','P','NOBS']
        for c in cols:
            col_name=mdf.loc[i,c]
            if pd.isna(col_name):
                if c=='snp':
                    df[c]=df[mdf.loc[i,"chr"]]+'_'+df[mdf.loc[i,"bp"]]
                if c=='n_eff':
                    df[c]=str(mdf.loc[i,"sample_size"]).strip()
            else:
                df[c]=df[col_name]
        df=df[cols]
        df.columns=final_cols
        def clean_chr(val):
            val = str(val).lower().replace('chr', '')
            if val == 'x':
                return '23'
            elif val == 'y':
                return '24'
            elif val in ['m', 'mt']:
                return '25'
            return val
        df['CHR'] = df['CHR'].apply(clean_chr)
        if pd.isna(mdf.loc[i,'snp']):
            df = pd.merge(df, snp_map, on=['CHR', 'BP'], how='left')
            df['SNP']=df['rsid']
            df.drop(columns=['rsid'], inplace=True)
        df.to_csv(out_path,sep='\t',index=False,lineterminator='\n',na_rep='.')
        log(f'{abbr}: save {df.shape[0]} snps')


def run_magma(nt=5):
    out_dir=f'{PROJ_RESULT_DIR}/magma'
    make_dir(out_dir)
    magma_bin = f'{LOCAL_DIR}/lib/magma/magma'
    cmds = []
    MAGMA_GWAS_dir=f'{PROJ_DATA_DIR}/gwas'
    meta_path=f'{MAGMA_GWAS_dir}/gwas.meta.xlsx'
    mdf=pd.read_excel(meta_path)
    for i in mdf.index:
        pheno=mdf.loc[i,'abbr']
        # if pheno not in ['PBC']:
        #     continue
        pop=str(mdf.loc[i,'pop']).strip().lower()
        magma_ref_geno=magma_ref_genos[pop]
        magma_gwas_path = f'{MAGMA_GWAS_dir}/{pheno}.magma.gwas.tsv'
        result_prefix = f'{out_dir}/{pheno}'
        make_dir(os.path.dirname(result_prefix))
        cmd1 = f'{magma_bin} --annotate --snp-loc {magma_gwas_path} --gene-loc {magma_gene_annot}' \
               f' --out {result_prefix}'
        cmd2 = f'{magma_bin} --bfile {magma_ref_geno} --gene-annot {result_prefix}.genes.annot ' \
               f'--pval {magma_gwas_path} ncol=NOBS --out {result_prefix}'
        cmd = f'{cmd1} && {cmd2}'
        cmds.append(cmd)
    batch_shell_plus(cmds, nt)
    pass

def make_scdrs_gs_file(top_n_genes=1000):
    magma_out_dir=f'{PROJ_RESULT_DIR}/magma'
    scdrs_gs_path=f'{PROJ_RESULT_DIR}/scdrs/gs/processed_geneset.gs'
    gdf = pd.read_table(magma_gene_annot, header=None)
    id_gene = {}
    for i in gdf.index:
        id_gene[gdf.loc[i, 0]] = gdf.loc[i, 5]
    datas=[]
    for f in os.listdir(magma_out_dir):
        if f.endswith('.genes.out'):
            pheno=f.split('.')[0]
            df=pd.read_table(f'{magma_out_dir}/{f}',sep='\s+')
            df.sort_values(by=['ZSTAT'],ascending=False,inplace=True,ignore_index=True)
            gs=[]
            for i in df.index[:top_n_genes]:
                gene=id_gene[df.loc[i,'GENE']]
                zs=df.loc[i,'ZSTAT']
                gs.append(f'{gene}:{zs}')
            datas.append([pheno,','.join(gs)])
    fdf=pd.DataFrame(datas,columns=['TRAIT','GENESET'])
    make_dir(os.path.dirname(scdrs_gs_path))
    fdf.to_csv(scdrs_gs_path,index=False,sep='\t',lineterminator='\n')


def main():
    make_magam_gwas_format()
    run_magma()
    make_scdrs_gs_file()

if __name__ == '__main__':
    main()
