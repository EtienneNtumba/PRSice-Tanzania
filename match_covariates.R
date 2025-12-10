# Read QC-passed samples
qc_samples <- read.table("data_without_qcfinal.qc.fam", header=FALSE)
colnames(qc_samples) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

# Create sample ID
qc_samples$ID <- paste(qc_samples$FID, qc_samples$IID)

# Read existing full covariates
full_cov <- read.table("data_without_qcfinal_covariates.txt", header=TRUE)
full_cov$ID <- paste(full_cov$FID, full_cov$IID)

# Extract only QC-passed samples
matched_cov <- full_cov[full_cov$ID %in% qc_samples$ID, ]
matched_cov$ID <- NULL

# Save
write.table(matched_cov, "data_without_qcfinal_qc_covariates.txt", 
            row.names=FALSE, quote=FALSE, sep=" ")

cat("Original covariates samples:", nrow(full_cov), "\n")
cat("QC-passed samples:", nrow(qc_samples), "\n")
cat("Matched samples:", nrow(matched_cov), "\n")

if(nrow(matched_cov) == 0) {
  cat("ERROR: No matching samples found!\n")
  cat("\nFirst 5 QC samples:\n")
  print(head(qc_samples[,1:2], 5))
  cat("\nFirst 5 full covariate samples:\n")
  print(head(full_cov[,1:2], 5))
}
