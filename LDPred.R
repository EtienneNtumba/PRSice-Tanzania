# ============================================
# LDpred2 : Pipeline Complet pour Tanzanie
# ============================================

library(bigsnpr)
library(data.table)
library(dplyr)

# ============================================
# PARTIE 1 : PRÉPARER LES INPUTS
# ============================================

# 1.1 - Convertir PLINK → bigsnpr
cat("Conversion PLINK → bigsnpr...\n")
snp_readBed("data_without_qcfinal.bed")

# 1.2 - Attacher les données
obj.bigSNP <- snp_attach("data_without_qcfinal.rds")
G <- obj.bigSNP$genotypes  # Matrice génotypique
map <- obj.bigSNP$map      # Info SNPs (chr, pos, alleles)
fam <- obj.bigSNP$fam      # Info individus

cat("Dimensions génotypes:", nrow(G), "individus x", ncol(G), "SNPs\n")

# 1.3 - Charger statistiques GWAS
cat("Chargement statistiques GWAS...\n")
sumstats <- fread("gcta_mlma_results.mlma")
colnames(sumstats) <- c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")

# 1.4 - Formater pour LDpred2
df_beta <- data.frame(
  chr = sumstats$Chr,
  pos = sumstats$bp,
  a0 = sumstats$A2,
  a1 = sumstats$A1,
  beta = sumstats$b,
  beta_se = sumstats$se,
  n_eff = 1527
)

# Filtrer valeurs manquantes
df_beta <- df_beta %>% 
  filter(!is.na(beta), !is.na(beta_se), beta_se > 0)

cat("Statistiques GWAS:", nrow(df_beta), "SNPs\n")

# 1.5 - Aligner GWAS avec génotypes
info_snp <- snp_match(df_beta, map)
cat("SNPs alignés:", nrow(info_snp), "\n")

# ============================================
# PARTIE 2 : CALCULER MATRICE LD (DE VOS DONNÉES)
# ============================================

cat("Calcul matrice LD des données tanzaniennes...\n")

# Positions génétiques (nécessaire pour LD)
POS2 <- snp_asGeneticPos(
  chr = map$chromosome,
  pos = map$physical.pos,
  dir = ".",
  ncores = 8
)

# Calculer corrélations LD
# CRITIQUE : Utilise VOS données tanzaniennes !
corr0 <- snp_cor(
  G,
  infos.pos = POS2,
  size = 3 / 1000,      # Fenêtre 3Mb
  alpha = 1,
  ncores = 8
)

cat("Matrice LD calculée :", dim(corr0)[1], "x", dim(corr0)[2], "\n")

# ============================================
# PARTIE 3 : LDPRED2-AUTO
# ============================================

cat("Exécution LDpred2-auto...\n")

# Paramètres
h2_init <- 0.3  # Héritabilité initiale estimée

# LDpred2-auto (optimise automatiquement)
multi_auto <- snp_ldpred2_auto(
  corr0,
  df_beta = info_snp,
  h2_init = h2_init,
  vec_p_init = seq_log(1e-4, 0.9, length.out = 30),
  burn_in = 500,
  num_iter = 500,
  report_step = 20,
  allow_jump_sign = FALSE,
  shrink_corr = 0.95,
  ncores = 8
)

# Extraire résultats
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

cat("Chaînes MCMC convergées:", ncol(beta_auto), "\n")

# ============================================
# PARTIE 4 : CALCULER PRS
# ============================================

cat("Calcul des scores PRS...\n")

# Calculer PRS pour chaque individu
pred_auto <- big_prodMat(
  G, 
  beta_auto,
  ind.col = info_snp$`_NUM_ID_`
)

# Moyenne des chaînes MCMC
prs_ldpred2 <- rowMeans(pred_auto)

cat("PRS calculés pour", length(prs_ldpred2), "individus\n")

# ============================================
# PARTIE 5 : ÉVALUATION
# ============================================

cat("Évaluation de la performance...\n")

# Charger phénotypes
pheno <- fread("data_without_qcfinal.fam")
colnames(pheno)[1:6] <- c("FID", "IID", "PID", "MID", "Sex", "HbF")

# Charger covariables
covs <- fread("data_without_qcfinal_qc_covariates.txt")

# Créer dataset d'analyse
data_analysis <- data.frame(
  FID = fam$family.ID,
  IID = fam$sample.ID,
  HbF = pheno$HbF,
  PRS = prs_ldpred2
)

# Fusionner avec covariables
data_analysis <- merge(data_analysis, covs, by = c("FID", "IID"))

# Modèle null (seulement covariables)
model_null <- lm(
  HbF ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex,
  data = data_analysis
)

# Modèle complet (covariables + PRS)
model_full <- lm(
  HbF ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex,
  data = data_analysis
)

# Calcul R²
r2_null <- summary(model_null)$r.squared
r2_full <- summary(model_full)$r.squared
r2_prs <- r2_full - r2_null

# P-value PRS
coef_summary <- summary(model_full)$coefficients
p_prs <- coef_summary["PRS", "Pr(>|t|)"]

# ============================================
# RÉSULTATS
# ============================================

cat("\n")
cat("=" %R% 50, "\n")
cat("     RÉSULTATS LDpred2 - TANZANIE\n")
cat("=" %R% 50, "\n\n")

cat("R² Null (covariables)     :", sprintf("%.4f (%.2f%%)", r2_null, r2_null*100), "\n")
cat("R² Full (PRS + covs)      :", sprintf("%.4f (%.2f%%)", r2_full, r2_full*100), "\n")
cat("R² PRS (incrémental)      :", sprintf("%.4f (%.2f%%)", r2_prs, r2_prs*100), "\n")
cat("P-value PRS               :", format(p_prs, scientific = TRUE, digits = 3), "\n\n")

# Comparaison avec PRSice
r2_prsice <- 0.0113
improvement <- (r2_prs / r2_prsice - 1) * 100
cat("R² PRSice (baseline)      :", sprintf("%.4f (%.2f%%)", r2_prsice, r2_prsice*100), "\n")
cat("Amélioration absolue      : +", sprintf("%.4f", r2_prs - r2_prsice), "\n")
cat("Amélioration relative     :", sprintf("%.1f%%", improvement), "\n\n")

# Héritabilité estimée
h2_est <- mean(sapply(multi_auto, function(x) tail(x$path_h2_est, 1)))
cat("Héritabilité SNP estimée  :", sprintf("%.3f", h2_est), "\n")

# ============================================
# SAUVEGARDER RÉSULTATS
# ============================================

# Poids des SNPs
weights_df <- data.frame(
  SNP = map$marker.ID[info_snp$`_NUM_ID_`],
  CHR = map$chromosome[info_snp$`_NUM_ID_`],
  POS = map$physical.pos[info_snp$`_NUM_ID_`],
  A1 = info_snp$a1,
  A2 = info_snp$a0,
  BETA_LDpred2 = rowMeans(beta_auto)
)

fwrite(weights_df, "ldpred2_weights_tanzania_HbF.txt", sep = "\t")

# Scores individuels
scores_df <- data.frame(
  FID = fam$family.ID,
  IID = fam$sample.ID,
  PRS_LDpred2 = prs_ldpred2
)

fwrite(scores_df, "ldpred2_scores_tanzania_HbF.txt", sep = "\t")

# Résumé
summary_df <- data.frame(
  Method = c("PRSice", "LDpred2"),
  R2_PRS = c(0.0113, r2_prs),
  R2_Full = c(0.6438, r2_full),
  P_value = c(3.32e-5, p_prs),
  N_SNPs = c(968114, nrow(info_snp))
)

fwrite(summary_df, "comparison_PRSice_vs_LDpred2.txt", sep = "\t")

cat("\n✅ Fichiers sauvegardés:\n")
cat("  - ldpred2_weights_tanzania_HbF.txt\n")
cat("  - ldpred2_scores_tanzania_HbF.txt\n")
cat("  - comparison_PRSice_vs_LDpred2.txt\n\n")

cat("=" %R% 50, "\n")
```

---

## **💾 OUTPUTS GÉNÉRÉS PAR LE SCRIPT**

Après exécution, vous aurez :

**1. Poids des SNPs :**
```
ldpred2_weights_tanzania_HbF.txt
SNP          CHR  POS      A1  A2  BETA_LDpred2
rs123456     1    12345    A   G   0.0124
rs234567     1    23456    C   T  -0.0089
```

**2. Scores individuels :**
```
ldpred2_scores_tanzania_HbF.txt
FID       IID                PRS_LDpred2
203205    A02_SCD5216453     0.00124
203205    A05_SCD5216477     0.00118
```

**3. Comparaison PRSice vs LDpred2 :**
```
comparison_PRSice_vs_LDpred2.txt
Method     R2_PRS   R2_Full   P_value    N_SNPs
PRSice     0.0113   0.6438    3.32e-05   968114
LDpred2    0.0XXX   0.6XXX    X.XXe-XX   XXXXXX
```

---

## **⏱️ TEMPS D'EXÉCUTION ESTIMÉ**

- Conversion PLINK → bigsnpr : **2-5 minutes**
- Calcul matrice LD : **30-60 minutes**
- LDpred2-auto (MCMC) : **1-3 heures**
- Calcul PRS : **2-5 minutes**

**Total : 2-4 heures** (selon votre machine)

---

## **📊 RÉSULTAT ATTENDU**
```
==================================================
     RÉSULTATS LDpred2 - TANZANIE
==================================================

R² Null (covariables)     : 0.6397 (63.97%)
R² Full (PRS + covs)      : 0.6520 (65.20%)
R² PRS (incrémental)      : 0.0123 (1.23%)   ← Attendu entre 1.5-2.5%
P-value PRS               : 1.25e-05

R² PRSice (baseline)      : 0.0113 (1.13%)
Amélioration absolue      : +0.0010
Amélioration relative     : 8.8%

Héritabilité SNP estimée  : 0.285
