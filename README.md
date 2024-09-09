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
process_vcf("path/to/file.vcf")
## list of files
process_VCFs(vector_of_vcfs)
## All the VCF files in a folder
process_folder("path/to/folder/")
```





