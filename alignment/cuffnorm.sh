#!/bin/bash
#SBATCH -J cuffnorm
#SBATCH -o cuffnorm.%j.out
#SBATCH -e cuffnorm.%j.err
#SBATCH -p cn-short
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000cns

export PATH=/lustre1/lch3000_pkuhpc/liuyt/liuyf/python/Python-2.7.11:$PATH
#export PYTHONPATH=/lustre1/lch3000_pkuhpc/liuyt/liuyf/python/Python-2.7.11
#export PATH=$HOME/.local/bin:$PATH
#export PYTHON_EGG_CACHE=
export PATH=/lustre1/lch3000_pkuhpc/liuyt/liuyf/app/cufflinks-2.2.1.Linux_x86_64:$PATH

cuffnorm -no-update-check -p 20 \
         --use-sample-sheet \
         -o cuffnorm_D200709\
         /lustre1/lch3000_pkuhpc/liuyf/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf  \
         ./sample_sheet.txt
