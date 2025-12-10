#!/bin/bash

# ============================================
# Script d'aide pour LDpred2 Analysis
# ============================================

# Couleurs pour output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║         LDpred2 Analysis - Tanzania HbF                    ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""

# ============================================
# EXEMPLE D'UTILISATION BASIQUE
# ============================================

echo -e "${GREEN}📝 UTILISATION BASIQUE:${NC}"
echo ""
echo "Rscript ldpred2_analysis.R \\"
echo "  --bed data_without_qcfinal \\"
echo "  --gwas gcta_mlma_results.mlma \\"
echo "  --covs data_without_qcfinal_qc_covariates.txt \\"
echo "  --out tanzania_hbf_ldpred2 \\"
echo "  --n-gwas 1527 \\"
echo "  --ncores 8"
echo ""

# ============================================
# EXEMPLE AVEC TOUS LES PARAMÈTRES
# ============================================

echo -e "${GREEN}📝 UTILISATION COMPLÈTE (tous paramètres):${NC}"
echo ""
echo "Rscript ldpred2_analysis.R \\"
echo "  --bed data_without_qcfinal \\"
echo "  --gwas gcta_mlma_results.mlma \\"
echo "  --covs data_without_qcfinal_qc_covariates.txt \\"
echo "  --out tanzania_hbf_ldpred2 \\"
echo "  --n-gwas 1527 \\"
echo "  --h2-init 0.3 \\"
echo "  --ncores 8 \\"
echo "  --ld-window 3 \\"
echo "  --burn-in 500 \\"
echo "  --num-iter 500 \\"
echo "  --r2-prsice 0.0113 \\"
echo "  --verbose"
echo ""

# ============================================
# AIDE COMPLÈTE
# ============================================

echo -e "${GREEN}❓ AIDE COMPLÈTE:${NC}"
echo ""
echo "Rscript ldpred2_analysis.R --help"
echo ""

# ============================================
# ARGUMENTS DISPONIBLES
# ============================================

echo -e "${GREEN}📋 ARGUMENTS DISPONIBLES:${NC}"
echo ""
echo -e "${YELLOW}Arguments REQUIS:${NC}"
echo "  -b, --bed FILE        Prefix fichiers PLINK (.bed/.bim/.fam)"
echo "  -g, --gwas FILE       Fichier statistiques GWAS (format GCTA)"
echo "  -c, --covs FILE       Fichier covariables (PCs + Sex)"
echo "  -n, --n-gwas INT      Taille échantillon GWAS"
echo ""

echo -e "${YELLOW}Arguments OPTIONNELS:${NC}"
echo "  -o, --out PREFIX      Prefix fichiers sortie [défaut: ldpred2_output]"
echo "  --h2-init FLOAT       Héritabilité initiale [défaut: 0.3]"
echo "  --ncores INT          Nombre de cores [défaut: 1]"
echo "  --ld-window FLOAT     Fenêtre LD en Mb [défaut: 3]"
echo "  --burn-in INT         Itérations burn-in [défaut: 500]"
echo "  --num-iter INT        Itérations MCMC [défaut: 500]"
echo "  --r2-prsice FLOAT     R² PRSice pour comparaison [défaut: 0.0113]"
echo "  --pheno-col INT       Colonne phénotype .fam [défaut: 6]"
echo "  -v, --verbose         Mode verbeux"
echo "  --no-plots            Ne pas générer graphiques"
echo ""

# ============================================
# FICHIERS DE SORTIE
# ============================================

echo -e "${GREEN}📂 FICHIERS DE SORTIE:${NC}"
echo ""
echo "  PREFIX_weights.txt            - Poids beta pour chaque SNP"
echo "  PREFIX_scores.txt             - Scores PRS individuels"
echo "  PREFIX_summary.txt            - Résumé comparatif"
echo "  PREFIX_results_detailed.txt   - Résultats détaillés"
echo "  PREFIX_convergence.pdf        - Graphiques convergence MCMC"
echo "  PREFIX_prs_distribution.pdf   - Distribution des scores"
echo ""

# ============================================
# VÉRIFICATIONS SYSTÈME
# ============================================

echo -e "${GREEN}🔍 VÉRIFICATIONS:${NC}"
echo ""

# Vérifier R
if command -v Rscript &> /dev/null; then
    R_VERSION=$(Rscript --version 2>&1 | head -n1)
    echo -e "  ✅ R installé: ${R_VERSION}"
else
    echo -e "  ${RED}❌ R non installé${NC}"
fi

# Vérifier packages R
echo ""
echo "Vérification des packages R..."
Rscript -e "
packages <- c('optparse', 'bigsnpr', 'data.table', 'dplyr')
installed <- packages %in% rownames(installed.packages())
for(i in 1:length(packages)) {
  if(installed[i]) {
    cat('  ✅', packages[i], '\n')
  } else {
    cat('  ❌', packages[i], '(NON INSTALLÉ)\n')
  }
}
" 2>/dev/null

echo ""

# ============================================
# INSTALLATION DES PACKAGES
# ============================================

echo -e "${GREEN}📦 INSTALLATION DES PACKAGES (si nécessaire):${NC}"
echo ""
echo "Rscript -e \""
echo "  install.packages('optparse')"
echo "  install.packages('data.table')"
echo "  install.packages('dplyr')"
echo "  install.packages('remotes')"
echo "  remotes::install_github('privefl/bigsnpr')"
echo "\""
echo ""

# ============================================
# TEMPS D'EXÉCUTION ESTIMÉ
# ============================================

echo -e "${GREEN}⏱️  TEMPS D'EXÉCUTION ESTIMÉ:${NC}"
echo ""
echo "  Conversion PLINK      : 2-5 minutes"
echo "  Calcul matrice LD     : 30-60 minutes"
echo "  LDpred2-auto (MCMC)   : 1-3 heures"
echo "  Calcul PRS            : 2-5 minutes"
echo "  ─────────────────────────────────"
echo "  TOTAL                 : 2-4 heures"
echo ""

# ============================================
# NOTES IMPORTANTES
# ============================================

echo -e "${GREEN}⚠️  NOTES IMPORTANTES:${NC}"
echo ""
echo "  1. Assurez-vous d'avoir au moins 16GB RAM"
echo "  2. Utilisez maximum 8 cores pour éviter surcharge"
echo "  3. Le script utilise le LD de VOS données (pas de panel externe)"
echo "  4. Formats GWAS supportés: GCTA MLMA (.mlma)"
echo "  5. Les fichiers PLINK doivent avoir extension .bed/.bim/.fam"
echo ""

# ============================================
# EXEMPLE COMPLET
# ============================================

echo -e "${GREEN}🚀 EXEMPLE D'EXÉCUTION COMPLET:${NC}"
echo ""
echo -e "${BLUE}# 1. Vérifier que les fichiers existent${NC}"
echo "ls data_without_qcfinal.bed"
echo "ls gcta_mlma_results.mlma"
echo "ls data_without_qcfinal_qc_covariates.txt"
echo ""

echo -e "${BLUE}# 2. Lancer l'analyse${NC}"
echo "Rscript ldpred2_analysis.R \\"
echo "  --bed data_without_qcfinal \\"
echo "  --gwas gcta_mlma_results.mlma \\"
echo "  --covs data_without_qcfinal_qc_covariates.txt \\"
echo "  --out tanzania_hbf_results \\"
echo "  --n-gwas 1527 \\"
echo "  --ncores 8 \\"
echo "  --verbose"
echo ""

echo -e "${BLUE}# 3. Vérifier les résultats${NC}"
echo "cat tanzania_hbf_results_summary.txt"
echo "cat tanzania_hbf_results_results_detailed.txt"
echo ""

# ============================================
# SUPPORT
# ============================================

echo -e "${GREEN}📧 SUPPORT:${NC}"
echo ""
echo "  Pour questions ou bugs, consultez la documentation"
echo "  GitHub: privefl/bigsnpr"
echo ""

echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
