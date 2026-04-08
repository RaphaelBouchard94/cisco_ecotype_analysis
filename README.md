# Genomic and phenotypic analysis of cisco (*Coregonus artedi*) ecotypes

This repository contains all custom scripts used for the analyses and figures presented in:

> Bouchard R. et al. (*in review*). [Indigenous knowledge reveals overlooked eco-evolutionary divergence in anadromous cisco]. *Nature Ecology & Evolution*.

---

## Study overview

This study investigates ecotype differentiation in cisco (*Coregonus artedi*) using a combination of low-coverage whole-genome sequencing (lcWGS), GTseq genotyping-by-sequencing, and phenotypic data (linear morphometrics, geometric morphometrics, life history traits). Analyses include population structure, genome scans for differentiation, linkage disequilibrium, detection of chromosomal inversions, population assignment, and candidate gene investigation.

---

## Repository structure

```
cisco_ecotype_analysis/
├── bash/                          # Bioinformatics pipeline (lcWGS)
│   ├── 01_config.sh               # Configuration file (paths, parameters)
│   ├── 02_list_bamfiles.sh        # Generate BAM file lists
│   ├── 03_saf_maf_gl_all_parallel_LL.sh  # SAF, MAF and genotype likelihoods (ANGSD)
│   ├── 03.1_make_sites_file.sh    # Create sites file
│   ├── 04_ngsParalog.sh           # Paralog detection (ngsParalog)
│   ├── 04.1_filter_based_on_lr.sh # Filter sites based on likelihood ratio
│   ├── 04.2_subset_canonical_sites.sh    # Subset to canonical sites
│   ├── 05.0_extract_singletons_parallel.sh  # Singleton extraction
│   ├── 05.1_merge.sh              # Merge singleton files
│   ├── 06_pca.sh                  # PCA (PCAngsd)
│   ├── 07_ngsadmix.sh             # Admixture analysis (NGSadmix)
│   ├── 08.0_generate_vcf.sh       # VCF generation
│   ├── 08.1_merger_vcf.sh         # Merge VCF files
│   ├── 08.2_run_ngsamova.sh       # AMOVA (ngsamova)
│   ├── 09_calculate_saf_per_chr.sh  # Per-chromosome SAF
│   ├── 09_calculate_maf_per_chr.sh  # Per-chromosome MAF
│   ├── 09_merge_saf_per_pop.sh    # Merge SAF by population
│   ├── 10_calculate_genome_wide_fst.sh   # Genome-wide Fst
│   ├── 10_print_by_pos_fst.sh     # Fst per position
│   ├── 12_folded_sfs.sh           # Folded site frequency spectrum
│   ├── 12_calculate_thetas.sh     # Theta statistics (ANGSD)
│   ├── 12_theta_per_window.sh     # Sliding window theta
│   ├── 13_dxy_chr34.sh            # Dxy on chromosome 34
│   ├── 13_dxy_per_window.sh       # Sliding window Dxy
│   ├── 14_ngsLD_parallel.sh       # Linkage disequilibrium (ngsLD)
│   ├── 14_window_ngsLD.sh         # Window-based LD
│   ├── 14_get_beagle_karyoaa.sh   # Beagle file for karyotype AA
│   ├── 14_get_beagle_karyobb.sh   # Beagle file for karyotype BB
│   ├── 14_get_beagle_notsummer.sh # Beagle file excluding summer samples
│   ├── 14_ngsLD_karyoaa.sh        # LD for karyotype AA
│   ├── 14_ngsLD_karyobb.sh        # LD for karyotype BB
│   ├── 14_ngsLD_notsummer.sh      # LD excluding summer samples
│   ├── 15_beagle_per_window.sh    # Beagle per window
│   ├── 15_beagle_per_window_imp.sh  # Imputed Beagle per window
│   ├── 15_generate_beagle_chr34_outlier.sh  # Beagle for chr34 outlier region
│   ├── 15_pca_outlier.sh          # PCA on outlier regions
│   ├── 15_sliding_pca.sh          # Sliding window PCA
│   ├── 16_het_chr34.sh            # Heterozygosity on chromosome 34
│   ├── 17_extract_fasta_sequence.sh  # FASTA extraction
│   ├── 18_gwas_length.sh          # GWAS for body length
│   ├── 19_depth.sh                # Sequencing depth calculation
│   └── [helper scripts]           # Python and R helper scripts used in pipeline
│
└── R/                             # Statistical analyses and figures
    ├── 00_map_rupert_bay.R        # Study area map
    ├── 01_linear_measurement_analysis.R     # Linear morphometric analysis
    ├── 02_geometric_morphometric_analysis.R # Geometric morphometric analysis
    ├── 02_geometric_morphometric_allometric_corr.R  # Allometric correction
    ├── 02_geometric_morphometric_arching_biais_corr.R  # Arching bias correction
    ├── 03_life_history_analysis.R # Age and growth analysis (von Bertalanffy)
    ├── 04_pca.R                   # Population genomic PCA
    ├── 04_admixture.R             # Admixture analysis visualization
    ├── 04_genome_scan_fst_clean.R # Genome scan for Fst outliers
    ├── 04_gwas_analysis_clean.R   # GWAS analysis and visualization
    ├── 04_gwas_pheno_file.R       # Phenotype file preparation for GWAS
    ├── 04_detect_outlier_windows_notsummer_rup.R  # Outlier window detection
    ├── 04_window_selection.R      # Outlier window selection
    ├── 05_inversion_analysis_pca.R  # Chromosomal inversion PCA
    ├── 05_pca_by_window.R         # Window-based PCA
    ├── 05_outlier_window_analysis.R  # Outlier window characterization
    ├── 05_het_chr36.R             # Heterozygosity analysis chromosome 36
    ├── 05_LD_Xplot.R              # LD cross-plot
    ├── 05_headtmap_ld_ch36.R      # LD heatmap chromosome 36
    ├── 05_gene_position.R         # Candidate gene position analysis
    ├── 05_mc4r_alignment.R        # MC4R protein alignment visualization
    ├── 05_genotype_length_analysis_gtseq.R  # GTseq genotype-length analysis
    ├── 06_genetic_assignment_analysis.R     # Population assignment (rubias)
    ├── 06_rubias_assignment_analysis.R      # Rubias assignment
    ├── 06_gtseq_genotyping.R      # GTseq genotyping pipeline
    ├── 06_gtseq_population_genetic.R  # GTseq population genetics
    ├── 06_low_cov_all_pop_pca.R   # All-population PCA (lcWGS)
    ├── 06_low_cov_all_pop_admixture.R  # All-population admixture (lcWGS)
    ├── 06_map_assignment_results.R  # Map of assignment results
    ├── 06_binomial_model_prob_kuukamek.R  # Binomial model for assignment probability
    ├── GTscore/                   # GTscore pipeline for GTseq data processing
    └── data/                      # Input data files (see Data availability)
```

## Usage

Scripts are numbered to reflect the order of execution. Start with `bash/01_config.sh` to set all paths and parameters for your system before running the pipeline. Each script contains a header comment describing its purpose and expected inputs/outputs.

The R scripts are largely independent and can be run in the order indicated by their numeric prefix. Scripts in `R/data/` contain the input data required by the R analyses.

---

## Data availability

Raw sequencing data (lcWGS and GTseq) are available on NCBI SRA under BioProject accession [XXXXX]. Intermediate files required to run the R scripts directly are available in `R/data/`.

---

## Citation

If you use these scripts, please cite:

> Bouchard R. et al. (*in review*). [Manuscript title]. *Nature Ecology & Evolution*.

---

## Contact

Raphaël Bouchard — raphael.bouchard.1@ulaval.ca  
