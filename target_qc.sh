#!/usr/bin/env bash

# Usage: ./target_qc.sh <target_prefix>
# Example: ./target_qc.sh my_target_dataset

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <target_prefix>"
    echo "  <target_prefix> = préfixe des fichiers PLINK binaires (.bed/.bim/.fam)"
    exit 1
fi

TARGET="$1"

echo "Running QC on target dataset: ${TARGET}"

plink \
  --bfile "${TARGET}" \
  --maf 0.05 \
  --mind 0.1 \
  --geno 0.1 \
  --hwe 1e-6 \
  --make-just-bim \
  --make-just-fam \
  --out "${TARGET}.qc"

echo
echo "QC finished."
echo "Generated files:"
echo "  ${TARGET}.qc.bim  -> SNPs passing QC"
echo "  ${TARGET}.qc.fam  -> individuals passing QC"
echo
echo "You can now add the following options to your PRSice command:"
echo "  --keep ${TARGET}.qc.fam --extract ${TARGET}.qc.bim"
echo
echo "For more details, see Marees et al. (2018) and the PRSice tutorial."

