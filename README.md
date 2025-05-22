# Variant Helpers

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
## Get standardised columns for Depth, RDP, ADP, AF, RDF, RDR, ADF, ADR
VCF <- process_counts(file, VCF, ...)
```

### Filtering

The resulting VCF files are just standard dataframes so you can filter them like normal. You can however use the filter_vars function to quickly apply some standard filters.

```r
## Test the default filters and see how it affects the data
filter_vars(VCF, plotStats = TRUE)
## Apply the default filters
filter_vars(VCF)
## set custom filters
filter_vars(VCF, vaf = 0.25, depth = 100, alt.depth = 10)
```

There are several filters that can be applied using this function

Argument | Default | Notes
--- | --- | ---
biotype | c("protein_coding") | [The type of feature the annotation is refering to](https://www.ensembl.org/info/genome/genebuild/biotypes.html)
impact | c("HIGH", "MODERATE") | [High and Moderate retain the non-synonymous variants](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
existing | TRUE | If TRUE then any variants that have an rs ID (from dbSNP) also need to have a COSMIC ID in the ID column
population | 0.01 | max frequency in 1000g and GnomAD
vaf | 0 | Tumour allele frequency
depth | 1 | minimum depth of coverage at position
alt.depth | 1 | minimum number of reads for the alternative allele

#### Filter preview

You can plot some summaries of variant filtering using the argument plotStats = TRUE, you can also return the plot objects or the underliying stats using returnPlots = TRUE or returnStats = TRUE

```r
VCF.filtered <- filter_vars(VCF, plotStats = TRUE)
plots <- filter_vars(VCF, returnPlots = TRUE)
stats <- filter_vars(VCF, returnStats = TRUE)
```

### Per caller plots

Prior to generating consensus calls you can take a look at the data across the callers.

You can see the number of calls per sample per caller using the command below

```r
## Plot the variants called by each vartiant caller.
variants_by_caller(VCF, category = Sample)
```

You can generate plots of the numerical columns (e.g. AF, DP, ADP) using the below

```r
## compare the AF profile across callers
stats_by_caller(TIER, category = Patient, Stat = AF)
```

This can let you check for caller specific issues with your variant calling.

### Consensus calling

When using multiple variant callers you can take the consensus of the calls from the different callers to filter out potential false calls.

#### Upset plots

Upset plots are a better way of comparing the overlap of multiple groups than Venn Diagrams. Venn Diagrams essentially stop being at all informative when you reach 5 + groups and are already pretty hard to quickly understand with 4 groups.

Upset graphs allow for a much better visualization of things shared between groups

```r
upset_by_caller(object = VCF, 
                category = Sample)
## You can specify the colours using a vector
upset_by_caller(object = VCF, 
                category = Sample,
                colours = c("vector", "of", "colours"))
## Even better - you can set up a named vector and use that throughout your code to standardise colours
sample_cols <- c("Sample1" = "#F48565", "Sample2" = "#85F465", "Sample3" = "#8565F4")
upset_by_caller(object = VCF, 
                category = Sample,
                colours = sample_cols)
```

#### Creating a consensus list of variants

Looking at the upset plot will give you a decent idea of how many callers to accept.

If you want to be really strict you can take only variants called by ALL of the callers you have used, or you could be very lax and use variants called by all tools. 

```r
## Add the consensus calls to the consensus slot
VCF <- add_consensus(VCF, min_caller = 1)
```
