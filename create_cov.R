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
