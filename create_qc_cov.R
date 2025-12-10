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
