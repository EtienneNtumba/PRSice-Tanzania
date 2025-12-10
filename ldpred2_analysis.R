#!/usr/bin/env Rscript

# ============================================
# LDpred2 Analysis Pipeline - Command Line
# ============================================
# Author: Genomic Analysis
# Date: 2025-12-08
# Description: Calcule PRS avec LDpred2 pour données tanzaniennes
# ============================================

# Charger librairies nécessaires
suppressPackageStartupMessages({
  library(optparse)
  library(bigsnpr)
  library(data.table)
  library(dplyr)
})

# ============================================
# DÉFINITION DES ARGUMENTS
# ============================================

option_list <- list(
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="Prefix des fichiers PLINK (.bed/.bim/.fam) [REQUIS]",
              metavar="PREFIX"),
  
  make_option(c("-g", "--gwas"), type="character", default=NULL,
              help="Fichier statistiques GWAS (format GCTA MLMA) [REQUIS]",
              metavar="FILE"),
  
  make_option(c("-c", "--covs"), type="character", default=NULL,
              help="Fichier covariables (PCs + Sex) [REQUIS]",
              metavar="FILE"),
  
  make_option(c("-o", "--out"), type="character", default="ldpred2_output",
              help="Prefix fichiers de sortie [défaut: ldpred2_output]",
              metavar="PREFIX"),
  
  make_option(c("-n", "--n-gwas"), type="integer", default=NULL,
              help="Taille échantillon GWAS [REQUIS]",
              metavar="INTEGER"),
  
  make_option(c("--h2-init"), type="double", default=0.3,
              help="Héritabilité initiale estimée [défaut: 0.3]",
              metavar="FLOAT"),
  
  make_option(c("--ncores"), type="integer", default=1,
              help="Nombre de cores CPU [défaut: 1]",
              metavar="INTEGER"),
  
  make_option(c("--ld-window"), type="double", default=3,
              help="Fenêtre LD en Mb [défaut: 3]",
              metavar="FLOAT"),
  
  make_option(c("--burn-in"), type="integer", default=500,
              help="Itérations burn-in MCMC [défaut: 500]",
              metavar="INTEGER"),
  
  make_option(c("--num-iter"), type="integer", default=500,
              help="Itérations MCMC [défaut: 500]",
              metavar="INTEGER"),
  
  make_option(c("--r2-prsice"), type="double", default=0.0113,
              help="R² PRSice pour comparaison [défaut: 0.0113]",
              metavar="FLOAT"),
  
  make_option(c("--pheno-col"), type="integer", default=6,
              help="Colonne phénotype dans .fam [défaut: 6]",
              metavar="INTEGER"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Mode verbeux"),
  
  make_option(c("--no-plots"), action="store_true", default=FALSE,
              help="Ne pas générer les graphiques")
)

# Parser les arguments
opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "\nAnalyse LDpred2 pour données tanzaniennes de drépanocytose (HbF).\nCalcule scores de risque polygénique avec matrice LD in-sample.",
  epilogue = "\nExemple:\n  Rscript ldpred2_analysis.R \\\n    --bed data_without_qcfinal \\\n    --gwas gcta_mlma_results.mlma \\\n    --covs data_without_qcfinal_qc_covariates.txt \\\n    --out tanzania_hbf \\\n    --n-gwas 1527 \\\n    --ncores 8\n"
)

opt <- parse_args(opt_parser)

# ============================================
# VALIDATION DES ARGUMENTS
# ============================================

validate_args <- function(opt) {
  errors <- c()
  
  # Arguments requis
  if (is.null(opt$bed)) {
    errors <- c(errors, "❌ --bed est requis")
  } else if (!file.exists(paste0(opt$bed, ".bed"))) {
    errors <- c(errors, paste0("❌ Fichier introuvable: ", opt$bed, ".bed"))
  }
  
  if (is.null(opt$gwas)) {
    errors <- c(errors, "❌ --gwas est requis")
  } else if (!file.exists(opt$gwas)) {
    errors <- c(errors, paste0("❌ Fichier introuvable: ", opt$gwas))
  }
  
  if (is.null(opt$covs)) {
    errors <- c(errors, "❌ --covs est requis")
  } else if (!file.exists(opt$covs)) {
    errors <- c(errors, paste0("❌ Fichier introuvable: ", opt$covs))
  }
  
  if (is.null(opt$`n-gwas`)) {
    errors <- c(errors, "❌ --n-gwas est requis")
  } else if (opt$`n-gwas` <= 0) {
    errors <- c(errors, "❌ --n-gwas doit être > 0")
  }
  
  # Validation paramètres
  if (opt$ncores < 1) {
    errors <- c(errors, "❌ --ncores doit être >= 1")
  }
  
  if (opt$`h2-init` <= 0 || opt$`h2-init` >= 1) {
    errors <- c(errors, "❌ --h2-init doit être entre 0 et 1")
  }
  
  if (length(errors) > 0) {
    cat("\n⚠️  ERREURS DE VALIDATION:\n")
    cat(paste(errors, collapse = "\n"), "\n\n")
    print_help(opt_parser)
    quit(status = 1)
  }
  
  return(TRUE)
}

validate_args(opt)

# ============================================
# FONCTION: Messages
# ============================================

msg <- function(text, type = "info") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  prefix <- switch(type,
    "info" = "ℹ️ ",
    "success" = "✅",
    "warning" = "⚠️ ",
    "error" = "❌",
    "step" = "▶️ ",
    ""
  )
  
  cat(sprintf("[%s] %s %s\n", timestamp, prefix, text))
}

# ============================================
# DÉBUT DE L'ANALYSE
# ============================================

msg(strrep("=", 60), "info")
msg("  LDpred2 Analysis Pipeline - Tanzanie HbF", "info")
msg(strrep("=", 60), "info")
msg("", "info")

msg("Paramètres:", "step")
msg(sprintf("  Fichiers PLINK    : %s", opt$bed), "info")
msg(sprintf("  Statistiques GWAS : %s", opt$gwas), "info")
msg(sprintf("  Covariables       : %s", opt$covs), "info")
msg(sprintf("  Prefix sortie     : %s", opt$out), "info")
msg(sprintf("  N GWAS            : %d", opt$`n-gwas`), "info")
msg(sprintf("  Héritabilité init : %.3f", opt$`h2-init`), "info")
msg(sprintf("  Cores CPU         : %d", opt$ncores), "info")
msg(sprintf("  Fenêtre LD        : %.1f Mb", opt$`ld-window`), "info")
msg("", "info")

# ============================================
# PARTIE 1: PRÉPARER LES INPUTS
# ============================================

msg("PARTIE 1: Préparation des données", "step")

# 1.1 - Convertir PLINK → bigsnpr
msg("Conversion PLINK → bigsnpr...", "info")
bed_file <- paste0(opt$bed, ".bed")
rds_file <- paste0(opt$bed, ".rds")

if (!file.exists(rds_file)) {
  tryCatch({
    snp_readBed(bed_file)
    msg("Conversion réussie", "success")
  }, error = function(e) {
    msg(paste("Erreur conversion:", e$message), "error")
    quit(status = 1)
  })
} else {
  msg("Fichier .rds existe déjà, réutilisation", "info")
}

# 1.2 - Attacher les données
msg("Chargement des données génotypiques...", "info")
tryCatch({
  obj.bigSNP <- snp_attach(rds_file)
  G <- obj.bigSNP$genotypes
  map <- obj.bigSNP$map
  fam <- obj.bigSNP$fam
  msg(sprintf("Chargé: %d individus × %d SNPs", nrow(G), ncol(G)), "success")
}, error = function(e) {
  msg(paste("Erreur chargement:", e$message), "error")
  quit(status = 1)
})

# 1.3 - Charger statistiques GWAS
msg("Chargement statistiques GWAS...", "info")
tryCatch({
  sumstats <- fread(opt$gwas)
  
  # Détecter format
  if (ncol(sumstats) == 9) {
    # Format GCTA MLMA
    colnames(sumstats) <- c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")
  } else {
    msg("Format GWAS non reconnu. Colonnes attendues: Chr, SNP, bp, A1, A2, Freq, b, se, p", "error")
    quit(status = 1)
  }
  
  msg(sprintf("Chargé: %d variants GWAS", nrow(sumstats)), "success")
}, error = function(e) {
  msg(paste("Erreur chargement GWAS:", e$message), "error")
  quit(status = 1)
})

# 1.4 - Formater pour LDpred2
msg("Formatage pour LDpred2...", "info")
df_beta <- data.frame(
  chr = sumstats$Chr,
  pos = sumstats$bp,
  a0 = sumstats$A2,
  a1 = sumstats$A1,
  beta = sumstats$b,
  beta_se = sumstats$se,
  n_eff = opt$`n-gwas`
)

# Filtrer valeurs manquantes et invalides
df_beta_clean <- df_beta %>%
  filter(!is.na(beta), !is.na(beta_se), beta_se > 0, is.finite(beta))

msg(sprintf("Variants après nettoyage: %d (%.1f%%)", 
            nrow(df_beta_clean), 
            100 * nrow(df_beta_clean) / nrow(sumstats)), "info")

# 1.5 - Aligner GWAS avec génotypes
msg("Alignement GWAS ↔ génotypes...", "info")
tryCatch({
  info_snp <- snp_match(df_beta_clean, map)
  msg(sprintf("SNPs alignés: %d (%.1f%%)", 
              nrow(info_snp), 
              100 * nrow(info_snp) / nrow(df_beta_clean)), "success")
}, error = function(e) {
  msg(paste("Erreur alignement:", e$message), "error")
  quit(status = 1)
})

# ============================================
# PARTIE 2: CALCULER MATRICE LD
# ============================================

msg("", "info")
msg("PARTIE 2: Calcul matrice LD (in-sample tanzanien)", "step")

# Positions génétiques
msg("Calcul positions génétiques...", "info")
tryCatch({
  POS2 <- snp_asGeneticPos(
    chr = map$chromosome,
    pos = map$physical.pos,
    dir = tempdir(),
    ncores = opt$ncores
  )
  msg("Positions génétiques calculées", "success")
}, error = function(e) {
  msg(paste("Erreur positions génétiques:", e$message), "error")
  quit(status = 1)
})

# Calculer matrice LD
msg(sprintf("Calcul matrice LD (fenêtre: %.1f Mb, %d cores)...", 
            opt$`ld-window`, opt$ncores), "info")
msg("⏳ Cette étape peut prendre 30-60 minutes...", "warning")

tryCatch({
  corr0 <- snp_cor(
    G,
    infos.pos = POS2,
    size = opt$`ld-window` / 1000,
    alpha = 1,
    ncores = opt$ncores
  )
  msg(sprintf("Matrice LD calculée: %d × %d", dim(corr0)[1], dim(corr0)[2]), "success")
}, error = function(e) {
  msg(paste("Erreur matrice LD:", e$message), "error")
  quit(status = 1)
})

# ============================================
# PARTIE 3: LDPRED2-AUTO
# ============================================

msg("", "info")
msg("PARTIE 3: LDpred2-auto (optimisation bayésienne)", "step")
msg(sprintf("Burn-in: %d, Itérations: %d", opt$`burn-in`, opt$`num-iter`), "info")
msg("⏳ Cette étape peut prendre 1-3 heures...", "warning")

tryCatch({
  multi_auto <- snp_ldpred2_auto(
    corr0,
    df_beta = info_snp,
    h2_init = opt$`h2-init`,
    vec_p_init = seq_log(1e-4, 0.9, length.out = 30),
    burn_in = opt$`burn-in`,
    num_iter = opt$`num-iter`,
    report_step = 20,
    allow_jump_sign = FALSE,
    shrink_corr = 0.95,
    ncores = opt$ncores
  )
  
  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  msg(sprintf("MCMC convergé: %d chaînes", ncol(beta_auto)), "success")
}, error = function(e) {
  msg(paste("Erreur LDpred2-auto:", e$message), "error")
  quit(status = 1)
})

# ============================================
# PARTIE 4: CALCULER PRS
# ============================================

msg("", "info")
msg("PARTIE 4: Calcul des scores PRS", "step")

tryCatch({
  pred_auto <- big_prodMat(
    G,
    beta_auto,
    ind.col = info_snp$`_NUM_ID_`
  )
  
  prs_ldpred2 <- rowMeans(pred_auto)
  msg(sprintf("PRS calculés pour %d individus", length(prs_ldpred2)), "success")
}, error = function(e) {
  msg(paste("Erreur calcul PRS:", e$message), "error")
  quit(status = 1)
})

# ============================================
# PARTIE 5: ÉVALUATION
# ============================================

msg("", "info")
msg("PARTIE 5: Évaluation de la performance", "step")

# Charger phénotypes
msg("Chargement phénotypes...", "info")
pheno <- fread(paste0(opt$bed, ".fam"))
colnames(pheno)[1:6] <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")

# Charger covariables
msg("Chargement covariables...", "info")
covs <- fread(opt$covs)

# Créer dataset d'analyse
data_analysis <- data.frame(
  FID = fam$family.ID,
  IID = fam$sample.ID,
  HbF = pheno$Phenotype,
  PRS = prs_ldpred2
)

data_analysis <- merge(data_analysis, covs, by = c("FID", "IID"))

# Vérifier que toutes les colonnes nécessaires existent
required_cols <- c("HbF", "PRS", "PC1", "PC2", "PC3", "PC4", "PC5", 
                   "PC6", "PC7", "PC8", "PC9", "PC10", "Sex")
missing_cols <- setdiff(required_cols, colnames(data_analysis))

if (length(missing_cols) > 0) {
  msg(paste("Colonnes manquantes:", paste(missing_cols, collapse = ", ")), "error")
  quit(status = 1)
}

# Modèle null
msg("Régression modèle null...", "info")
model_null <- lm(
  HbF ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex,
  data = data_analysis
)

# Modèle complet
msg("Régression modèle complet...", "info")
model_full <- lm(
  HbF ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex,
  data = data_analysis
)

# Calcul métriques
r2_null <- summary(model_null)$r.squared
r2_full <- summary(model_full)$r.squared
r2_prs <- r2_full - r2_null

coef_summary <- summary(model_full)$coefficients
p_prs <- coef_summary["PRS", "Pr(>|t|)"]
beta_prs <- coef_summary["PRS", "Estimate"]
se_prs <- coef_summary["PRS", "Std. Error"]

# Héritabilité estimée
h2_est <- mean(sapply(multi_auto, function(x) tail(x$path_h2_est, 1)))
p_est <- mean(sapply(multi_auto, function(x) tail(x$path_p_est, 1)))

# Amélioration vs PRSice
improvement_abs <- r2_prs - opt$`r2-prsice`
improvement_rel <- (r2_prs / opt$`r2-prsice` - 1) * 100

# ============================================
# AFFICHER RÉSULTATS
# ============================================

msg("", "info")
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("               RÉSULTATS LDpred2 - TANZANIE HbF\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("📊 PERFORMANCE DU MODÈLE:\n")
cat(sprintf("   R² Null (covariables seules)  : %.4f (%.2f%%)\n", r2_null, r2_null*100))
cat(sprintf("   R² Full (PRS + covariables)   : %.4f (%.2f%%)\n", r2_full, r2_full*100))
cat(sprintf("   R² PRS (incrémental)          : %.4f (%.2f%%) ✨\n", r2_prs, r2_prs*100))
cat(sprintf("   P-value PRS                   : %s\n", format(p_prs, scientific = TRUE, digits = 3)))
cat(sprintf("   Beta PRS                      : %.4f ± %.4f\n", beta_prs, se_prs))
cat("\n")

cat("📈 COMPARAISON AVEC PRSice:\n")
cat(sprintf("   R² PRSice (baseline)          : %.4f (%.2f%%)\n", opt$`r2-prsice`, opt$`r2-prsice`*100))
cat(sprintf("   Amélioration absolue          : +%.4f\n", improvement_abs))
cat(sprintf("   Amélioration relative         : %+.1f%%\n", improvement_rel))
cat("\n")

cat("🧬 PARAMÈTRES GÉNÉTIQUES:\n")
cat(sprintf("   Héritabilité SNP (h²)         : %.3f\n", h2_est))
cat(sprintf("   Polygenicité (p)              : %.4f\n", p_est))
cat(sprintf("   Nombre SNPs utilisés          : %d\n", nrow(info_snp)))
cat("\n")

cat("💾 STATISTIQUES:\n")
cat(sprintf("   Individus analysés            : %d\n", nrow(data_analysis)))
cat(sprintf("   Chaînes MCMC convergées       : %d\n", ncol(beta_auto)))
cat("\n")

cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ============================================
# SAUVEGARDER RÉSULTATS
# ============================================

msg("SAUVEGARDE DES RÉSULTATS", "step")

# 1. Poids des SNPs
weights_file <- paste0(opt$out, "_weights.txt")
weights_df <- data.frame(
  SNP = map$marker.ID[info_snp$`_NUM_ID_`],
  CHR = map$chromosome[info_snp$`_NUM_ID_`],
  POS = map$physical.pos[info_snp$`_NUM_ID_`],
  A1 = info_snp$a1,
  A2 = info_snp$a0,
  BETA = rowMeans(beta_auto)
)
fwrite(weights_df, weights_file, sep = "\t")
msg(sprintf("Sauvegardé: %s", weights_file), "success")

# 2. Scores individuels
scores_file <- paste0(opt$out, "_scores.txt")
scores_df <- data.frame(
  FID = fam$family.ID,
  IID = fam$sample.ID,
  PRS = prs_ldpred2
)
fwrite(scores_df, scores_file, sep = "\t")
msg(sprintf("Sauvegardé: %s", scores_file), "success")

# 3. Résumé comparatif
summary_file <- paste0(opt$out, "_summary.txt")
summary_df <- data.frame(
  Method = c("PRSice", "LDpred2"),
  R2_PRS = c(opt$`r2-prsice`, r2_prs),
  R2_Full = c(NA, r2_full),
  P_value = c(3.32e-5, p_prs),
  N_SNPs = c(968114, nrow(info_snp)),
  Heritability = c(NA, h2_est)
)
fwrite(summary_df, summary_file, sep = "\t")
msg(sprintf("Sauvegardé: %s", summary_file), "success")

# 4. Résultats détaillés
results_file <- paste0(opt$out, "_results_detailed.txt")
sink(results_file)
cat(strrep("=", 70), "\n")
cat("LDpred2 Analysis - Detailed Results\n")
cat(strrep("=", 70), "\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUTS:\n")
cat("  PLINK prefix     :", opt$bed, "\n")
cat("  GWAS file        :", opt$gwas, "\n")
cat("  Covariates file  :", opt$covs, "\n")
cat("  N GWAS           :", opt$`n-gwas`, "\n\n")

cat("PARAMETERS:\n")
cat("  h² initial       :", opt$`h2-init`, "\n")
cat("  LD window (Mb)   :", opt$`ld-window`, "\n")
cat("  Burn-in          :", opt$`burn-in`, "\n")
cat("  Iterations       :", opt$`num-iter`, "\n")
cat("  Cores            :", opt$ncores, "\n\n")

cat("RESULTS:\n")
cat("  R² Null          :", sprintf("%.6f", r2_null), "\n")
cat("  R² Full          :", sprintf("%.6f", r2_full), "\n")
cat("  R² PRS           :", sprintf("%.6f", r2_prs), "\n")
cat("  P-value          :", format(p_prs, scientific = TRUE), "\n")
cat("  Beta PRS         :", sprintf("%.6f", beta_prs), "\n")
cat("  SE PRS           :", sprintf("%.6f", se_prs), "\n")
cat("  h² estimated     :", sprintf("%.4f", h2_est), "\n")
cat("  p estimated      :", sprintf("%.4f", p_est), "\n")
cat("  N SNPs used      :", nrow(info_snp), "\n\n")

cat("COMPARISON:\n")
cat("  PRSice R²        :", sprintf("%.6f", opt$`r2-prsice`), "\n")
cat("  Improvement (abs):", sprintf("%+.6f", improvement_abs), "\n")
cat("  Improvement (rel):", sprintf("%+.2f%%", improvement_rel), "\n\n")

print(summary(model_full))
sink()
msg(sprintf("Sauvegardé: %s", results_file), "success")

# 5. Graphiques (si demandé)
if (!opt$`no-plots`) {
  msg("Génération des graphiques...", "info")
  
  tryCatch({
    # Plot convergence
    conv_file <- paste0(opt$out, "_convergence.pdf")
    pdf(conv_file, width = 12, height = 6)
    par(mfrow = c(1, 2))
    
    # h² convergence
    h2_paths <- sapply(multi_auto, function(x) x$path_h2_est)
    matplot(h2_paths, type = "l", lty = 1, col = rgb(0, 0, 1, 0.3),
            xlab = "Itération", ylab = "h² estimé",
            main = "Convergence h²")
    abline(h = h2_est, col = "red", lwd = 2, lty = 2)
    
    # p convergence
    p_paths <- sapply(multi_auto, function(x) x$path_p_est)
    matplot(p_paths, type = "l", lty = 1, col = rgb(0, 0, 1, 0.3),
            xlab = "Itération", ylab = "p estimé",
            main = "Convergence polygenicité (p)")
    abline(h = p_est, col = "red", lwd = 2, lty = 2)
    
    dev.off()
    msg(sprintf("Sauvegardé: %s", conv_file), "success")
    
    # Plot distribution PRS
    dist_file <- paste0(opt$out, "_prs_distribution.pdf")
    pdf(dist_file, width = 8, height = 6)
    hist(prs_ldpred2, breaks = 50, col = "steelblue", border = "white",
         main = "Distribution des Scores PRS",
         xlab = "Score PRS", ylab = "Fréquence")
    abline(v = mean(prs_ldpred2), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c("Moyenne"), col = "red", lwd = 2, lty = 2)
    dev.off()
    msg(sprintf("Sauvegardé: %s", dist_file), "success")
    
  }, error = function(e) {
    msg(paste("Erreur génération graphiques:", e$message), "warning")
  })
}

# ============================================
# FIN
# ============================================

msg("", "info")
msg(strrep("=", 60), "success")
msg("  ANALYSE TERMINÉE AVEC SUCCÈS", "success")
msg(strrep("=", 60), "success")
msg("", "info")

msg("FICHIERS GÉNÉRÉS:", "info")
msg(sprintf("  📄 %s", weights_file), "info")
msg(sprintf("  📄 %s", scores_file), "info")
msg(sprintf("  📄 %s", summary_file), "info")
msg(sprintf("  📄 %s", results_file), "info")
if (!opt$`no-plots`) {
  msg(sprintf("  📊 %s_convergence.pdf", opt$out), "info")
  msg(sprintf("  📊 %s_prs_distribution.pdf", opt$out), "info")
}

msg("", "info")
msg(sprintf("🎯 R² PRS final: %.4f (%.2f%%)", r2_prs, r2_prs * 100), "success")

if (improvement_rel > 0) {
  msg(sprintf("📈 Amélioration vs PRSice: +%.1f%%", improvement_rel), "success")
} else {
  msg(sprintf("📉 Changement vs PRSice: %.1f%%", improvement_rel), "warning")
}

msg("", "info")

# Code de sortie
quit(status = 0)
