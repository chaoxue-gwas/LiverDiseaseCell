# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/04/07 20:17
# @Function: parameters

import platform

# Common
DATA_DIR = '/home/xc/local/data'
if platform.system() == 'Windows':
    DATA_DIR = r'E:\WorkData\syncHPC\home\data'

PROJ_NAME='LiverDisCell'
Result_ver='20250407'
PROJ_DATA_DIR=f'{DATA_DIR}/project_data/{PROJ_NAME}'
PROJ_RESULT_DIR=f'{DATA_DIR}/projects/{PROJ_NAME}/{Result_ver}'

RAW_DATA_DIR=r'E:\ShareData\MiuShare\miu\projects\6.liver diseaseGWAS_scRNAseq\2.resourse'

## MAGMA resources
magma_resource_dir=f'{DATA_DIR}/resources/ToolResource/MAGMA'
magma_gene_annot = f'{magma_resource_dir}/NCBI37.3/NCBI37.3.gene.loc'
magma_ref_genos = {pop:f'{magma_resource_dir}/g1000_{pop}/g1000_{pop}' for pop in ['eur','eas']}

## variant db
snp_rsid_db_hg19=r'E:\WorkData\syncHPC\home\data\resources\VariantAnnotation\rsid\hg19.common.vcf.gz'