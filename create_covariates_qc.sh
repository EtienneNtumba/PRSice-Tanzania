#!/bin/bash

# Create covariates for QC-passed samples ONLY
WORKDIR=/mnt/lustre/groups/CBBI0818/TZ/data/data_target_prisce
PREFIX=data_without_qcfinal

cd $WORKDIR

echo "Creating covariates for QC-passed samples only..."

# Step 1: LD pruning using ONLY QC-passed samples
plink --bfile $PREFIX \
      --keep ${PREFIX}.qc.fam \
      --indep-pairwise 200 50 0.25 \
      --out ${PREFIX}_qc_pruned

# Step 2: PCA using ONLY QC-passed samples
plink --bfile $PREFIX \
      --keep ${PREFIX}.qc.fam \
      --extract ${PREFIX}_qc_pruned.prune.in \
      --pca 10 \
      --out ${PREFIX}_qc_pca

# Step 3: Create covariates file
cat > create_qc_cov.R << 'REOF'
# Read PCA results
pca <- read.table("data_without_qcfinal_qc_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Read QC FAM file for sex
fam <- read.table("data_without_qcfinal.qc.fam", header=FALSE)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

# Merge
covariates <- merge(pca, fam[,c("FID","IID","SEX")], by=c("FID","IID"))

# Reorder columns
covariates <- covariates[,c("FID", "IID", paste0("PC", 1:10), "SEX")]
colnames(covariates)[ncol(covariates)] <- "Sex"

# Save
write.table(covariates, "data_without_qcfinal_qc_covariates.txt", 
            row.names=FALSE, quote=FALSE, sep=" ")

cat("Covariates created: data_without_qcfinal_qc_covariates.txt\n")
cat("Samples:", nrow(covariates), "\n")

# Verify match with QC fam
qc_fam <- read.table("data_without_qcfinal.qc.fam")
cat("QC fam samples:", nrow(qc_fam), "\n")
cat("Match:", all(paste(covariates$FID, covariates$IID) %in% paste(qc_fam$V1, qc_fam$V2)), "\n")
REOF

Rscript create_qc_cov.R

echo "Done!"
head data_without_qcfinal_qc_covariates.txt
