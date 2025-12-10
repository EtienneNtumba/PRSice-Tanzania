#!/bin/bash

# Script to create covariates for data_without_qcfinal
# Generates: PC1-PC10 and Sex

WORKDIR=/mnt/lustre/groups/CBBI0818/TZ/data/data_target_prisce
PREFIX=data_without_qcfinal

cd $WORKDIR

echo "Step 1: LD pruning for PCA..."
plink --bfile $PREFIX \
      --keep ${PREFIX}.qc.fam \
      --extract ${PREFIX}.qc.bim \
      --indep-pairwise 200 50 0.25 \
      --out ${PREFIX}_pruned

echo "Step 2: Running PCA..."
plink --bfile $PREFIX \
      --keep ${PREFIX}.qc.fam \
      --extract ${PREFIX}_pruned.prune.in \
      --pca 10 \
      --out ${PREFIX}_pca

echo "Step 3: Creating covariates file..."
# Create R script to combine PCA and sex
cat > create_cov.R << 'REOF'
# Read PCA results
pca <- read.table("data_without_qcfinal_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Read FAM file for sex
fam <- read.table("data_without_qcfinal.fam", header=FALSE)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

# Merge
covariates <- merge(pca, fam[,c("FID","IID","SEX")], by=c("FID","IID"))

# Reorder columns: FID IID PC1-PC10 Sex
covariates <- covariates[,c("FID", "IID", paste0("PC", 1:10), "SEX")]
colnames(covariates)[ncol(covariates)] <- "Sex"

# Save
write.table(covariates, "data_without_qcfinal_covariates.txt", 
            row.names=FALSE, quote=FALSE, sep=" ")

cat("Covariates file created: data_without_qcfinal_covariates.txt\n")
cat("Number of samples:", nrow(covariates), "\n")
cat("Columns:", colnames(covariates), "\n")
REOF

Rscript create_cov.R

echo "Done! Covariates file created: ${PREFIX}_covariates.txt"
head ${PREFIX}_covariates.txt
