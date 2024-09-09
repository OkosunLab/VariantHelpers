# Variant Helpers

An r package containing some functions to assist in the analysis of variant call data.

## Quickstart

### Installation

```r
remotes::install.github("OkosunLab/VatiantHelpers")
```

### Process vcf files

Run the full variant processing pipeline

```r
## single file
process_vcf("path/to/file.vcf", normalPattern = "normal")
## list of files
process_VCFs(vector_of_vcfs, normalPattern = "normal")
## All the VCF files in a folder
process_folder("path/to/folder/", normalPattern = "normal")
```

Not all paired variant callers will name the normal file NORMAL, some use the file name. So the normalPattern option can be provided to find the normal column using the pattern and change the name to NORMAL this will mean the columns will be lined up if you bind several together (default for process_VCFs and process_folder)

## Overview

This package has a bunch of functions to help speed up analysis.

### VCF loading and processing

The basic processing pipeline will run the following steps:

```r
## Load file skipping the ## comment lines
VCF <- read_VCF(file)
## Split VEP annotations
VCF <- split_vep(VCF)
## Split the tumour column and name using the format column
VCF <- split_format(VCF, ...)
## Get standardised columns for Depth, RDP, ADP, AF
VCF <- process_counts(file, VCF, ...)
```




