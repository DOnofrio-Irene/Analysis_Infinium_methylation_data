# DNA Methylation Analysis Pipeline for Infinium Data

This pipeline is designed to analyze DNA methylation data generated from the Illumina HumanMethylation450 Beadchip for the DNA/RNA Dynamics  project (MSc Bioinformatics, University of Bologna).
DNA methylation has been identified to be widely associated to complex diseases. Among biological platforms to profile DNA methylation in human, the Illumina Infinium HumanMethylation450 BeadChip (450K) has been accepted as one of the most efficient technologies. 


## Installation and Dependencies
This pipeline was runned in R (4.3.0 (2023-04-21 ucrt)).
You will need to have the following software and packages installed:

 ```r
# BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")

# minfi
BiocManager::install("minfi")

# Illumina manifest
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# future.apply 
install.packages("future.apply")

# factoextra
install.packages("factoextra")

# qqman
install.packages("qqman")

# gplots
install.packages("gplots")
```
