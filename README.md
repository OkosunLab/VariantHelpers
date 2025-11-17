# Variant Helpers

See the [wiki](https://github.com/OkosunLab/VariantHelpers/wiki) for full documentation.

## Contents

1. [Introduction](#introduction)
2. [Quickstart](#quickstart)
    1. [Installation](#installation)
    2. [Process VCF files](#process-vcf-files)
3. [Overview](#overview)
    1. [VCF loading and processing](#vcf-loading-and-processing)
    2. [Filtering](#filtering)
    3. [Consensus calling](#Consensus-calling)

## Introduction

An r package containing some functions to assist in the analysis of variant call data.

***IMPORTANT NOTE***

I have been getting an issue with running this packages recently through the ondemand implementation of R studio on Apocrita getting the following error:

```r
Error in base::suppressWarnings(base::try("_", silent = TRUE))
  invalid use of pipe placeholder (<input>:1:0)
```

if you get this try turning of R diagnostics - this seems to be an issue with this: https://github.com/rstudio/rstudio/issues/14713

## Quickstart

### Installation

To install the package you can use the following:

```r
remotes::install.github("OkosunLab/VatiantHelpers")
```

However the package is currently private so you need to provide github some credentials to be able to download the package. The simplest way to do this is to generate a ghg keey for your account and pass it in the install call:

```r
remotes::install_github("OkosunLab/VariantHelpers",
                            auth_token = "ghg_yourkeygoeshere")
```

### Process VCF files

Run the full variant processing pipeline

```r
library(VariantHelpers)
## single file
vcf <- process_vcf("path/to/file.vcf", tumourPattern = "tumour", normalPattern = "normal")
## list of files
vcfs <- process_VCFs(vector_of_vcfs, tumourPattern = "tumour", normalPattern = "normal")
## All the VCF files in a folder
vcfs <- process_folder("path/to/folder/", tumourPattern = "tumour", normalPattern = "normal")
```

#### VariantHelpers class

The variant helpers class can contain the data from several different variant callers and the package contains several functions for use with this class.

```r
## Take a vector of different folders and load variants into the variant helper class
Paths <- list("Mutect2" = "path/to/Mutect2/", "Vardict" = "path/to/Vardict/")
VCF <- variant_helper_from_folder_paths(paths = Paths, tumourPattern = "tumour", normalPattern = "normal")
```

Not all paired variant callers will name the normal file NORMAL, some use the file name. So the normalPattern option can be provided to find the normal column using the pattern and change the name to NORMAL this will mean the columns will be lined up if you bind several together (default for process_VCFs and process_folder)
