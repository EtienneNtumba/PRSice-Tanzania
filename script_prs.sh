#!/bin/bash

Rscript PRSice.R \
    --dir . \
    --prsice ./PRSice_linux \
    --base /mnt/lustre/groups/CBBI0818/TZ/data/qc_results/gcta/gcta_mlma_results.mlma \
    --target /mnt/lustre/groups/CBBI0818/TZ/data/data_target_prisce/data_without_qcfinal \
    --keep /mnt/lustre/groups/CBBI0818/TZ/data/data_target_prisce/data_without_qcfinal.qc.fam \
    --extract /mnt/lustre/groups/CBBI0818/TZ/data/data_target_prisce/data_without_qcfinal.qc.bim \
    --thread 8 \
    --stat b \
    --beta \
    --binary-target F \
    --snp SNP \
    --chr Chr \
    --bp bp \
    --A1 A1 \
    --A2 A2 \
    --pvalue p \
    --cov /mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step5_final/covariates.txt \
    --cov-col PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Sex \
    --clump-kb 250 \
    --clump-r2 0.1 \
    --clump-p 1 \
    --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --all-score \
    --quantile 10 \
    --out data_without_qcfinal_PRS
