# DNA Methylation Analysis Pipeline for Infinium Data

This is the DNA/RNA Dynamics exam project (MSc Bioinformatics, University of Bologna), which aims to analyze DNA methylation data generated from the Illumina HumanMethylation450 Beadchip.
DNA methylation has been identified to be widely associated to complex diseases. Among biological platforms to profile DNA methylation in human, the Illumina Infinium HumanMethylation450 BeadChip (450K) has been accepted as one of the most efficient technologies. 

*DISCLAIMER*: *The choices in some passages (i.e. choice of the statistical test) were made for educational purposes only*


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


### Infinium HumanMethylation450K Manifest file
In some steps it is needed to check the Manifest file, which can be found on the [Illumina website](http://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html)

The manifest was then cleaned, by removing the control and the rs probes in the followinf way

```r
Illumina450Manifest_clean <- Illumina450Manifest[!Illumina450Manifest$IlmnID %in% notMappedToCHR$IlmnID,]
```

## WORKFLOW
