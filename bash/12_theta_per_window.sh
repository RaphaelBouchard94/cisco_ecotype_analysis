#!/bin/bash

OUT_DIR="12.1_theta_per_window"

# Merge notsummer population windows
echo "Merging notsummer windows..."
cat ${OUT_DIR}/chr36_outlier_notsummer_*.theta.thetas.idx.pestPG | \
  awk 'NR==1 || !/^#/' > ${OUT_DIR}/chr36_outlier_notsummer_merged.pestPG

# Merge rup population windows
echo "Merging rup windows..."
cat ${OUT_DIR}/chr36_outlier_rup_*.theta.thetas.idx.pestPG | \
  awk 'NR==1 || !/^#/' > ${OUT_DIR}/chr36_outlier_rup_merged.pestPG

echo "Merged files created!"
