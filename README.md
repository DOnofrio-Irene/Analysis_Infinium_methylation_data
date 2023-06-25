# DNA Methylation Analysis Pipeline for Infinium Data

This is the DNA/RNA Dynamics exam project (MSc Bioinformatics, University of Bologna), which aims to analyze DNA methylation data generated from the Illumina HumanMethylation450 Beadchip.
DNA methylation has been identified to be widely associated to complex diseases. Among biological platforms to profile DNA methylation in human, the Illumina Infinium HumanMethylation450 BeadChip (450K) has been accepted as one of the most efficient technologies. 

*DISCLAIMER*: *The choices in some passages were made for educational purposes only*

## WORKFLOW
* **Importing Data**
* **Quality Control**: Perform data quality assessment to identify any potential issues or biases in the raw data.
* **Between-array normalization**: apply normalization techniques to minimize technical variations and batch effects within and across samples.
* **Differential methylation analysis**: Identify differentially methylated regions (DMRs) or individual CpG sites associated with various conditions or phenotypes of interest.



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

# viridis
install.packages("viridis")
```
