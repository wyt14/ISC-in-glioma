# ISC-in-glioma
Transcriptome profiling reveals an immunosuppressive class in glioma with prognostic and immunotherapeutic relevance

## Overview
This project stores the primary raw data and analysis code used in this article. All data are obtained from public sources, specifically TCGA or CGGA, and are integrated for downstream analysis.

## Detailed file descriptions
```
├── 1.LGG                     	# Transcriptome expression data for LGG from the TCGA project  
│   ├── clinical.tsv          	# Clinical phenotype data for LGG  
│   └── rank.5.Rdata          	# Data after classifying LGG expression profiles into 5 types using NMF analysis  
├── 2.LGG.methylation         	# Methylation data for LGG, aligned with transcriptome data from TCGA  
│   ├── DMP_DE.Rdata                  # Differential methylation results  
│   ├── DMP.Rdata                     # Differential analysis results of methylation  
│   ├── LGG.methylation.group.Rdata   # Raw data for hierarchical clustering of methylation  
│ 
├── 3.GBM                     # Transcriptome expression data for GBM from the TCGA project  
│   ├── clinical.tsv          # Clinical phenotype data for GBM  
│   └── rank.5.Rdata          # Data after classifying GBM expression profiles into 5 types using NMF analysis  
├── 4.GBM.methylation         # Methylation data for GBM, aligned with transcriptome data from TCGA  
│   ├── DMP_DE.Rdata                  # Differential methylation results  
│   ├── DMP.Rdata                     # Differential analysis results of methylation  
│   ├── GBM.methylation.group.Rdata   # Raw data for hierarchical clustering of methylation  
│
├── 5.LGG.CNV           # Raw genotypic data for LGG, including different mutation types  
│   ├── All.LGG.maf     # Raw genotype data for LGG  
│   ├── LGG.CNV.rda     # CNV data obtained via GISTIC analysis for LGG  
│   └── LGG.SNP.rda     # SNP data obtained via GISTIC analysis for LGG  
├── 6.GBM.CNV           # Raw genotypic data for GBM, including different mutation types  
│   ├── All.GBM.maf     # Raw genotype data for GBM  
│   ├── GBM.CNV.rda     # CNV data obtained via GISTIC analysis for GBM  
│   └── GBM.SNP.rda     # SNP data obtained via GISTIC analysis for GBM
│
├── plot.GBM.R          # Visualization of analysis results for GBM transcriptome and methylation data  
├── plot.LGG.R          # Visualization of analysis results for LGG transcriptome and methylation data  
└── README.md
```

## Copyright and license
Copyright (c) 2025 Yuting Wang.
The files are licensed under TMU.
