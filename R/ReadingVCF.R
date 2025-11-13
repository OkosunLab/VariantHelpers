#' @description
#' A function fully process the VCF through the pipeline.
#'
#' @title process_vcf
#' @param file The path to the file to be read
#' @return A dataframe of the VCF skipping the commented lines
#' @keywords VCF
#' @importFrom tibble remove_rownames
#' @export
#' @examples
#'
#' process_vcf(file)

process_vcf <- function(file, ...) {
    ## Load file
    VCF <- read_VCF(file, ...)
    ## Split VEP annotations
    VCF <- split_vep(VCF, ...)
    ## Split format
    VCF <- split_format(VCF, ...)
    ## Get counts setup
    VCF <- process_counts(file, VCF, ...)
}

#' @description
#' A function to process a vector of VCF files
#'
#' @title process_VCFs
#' @param VCFs a vector of paths to VCF files
#' @param JoinAll bool of whether to join all the files together (default TRUE)
#' @return A dataframe of the VCF skipping the commented lines
#' @keywords VCF
#' @importFrom tibble remove_rownames
#' @export
#' @examples
#'
#' process_VCFs(file)

process_VCFs <- function(VCFs, JoinAll = TRUE, ...) {
    counter = 0
    print(paste0("Reading variant calls from ", length(VCFs), " files"))
    pb = txtProgressBar(min = 0,
                        max = length(VCFs),
                        initial = 0,
                        style = 3)
    VCFs <-
        lapply(VCFs, function(file) {
            counter <<- counter + 1
            vcf <- process_vcf(file, ...)
            setTxtProgressBar(pb,counter)
            vcf
        })
    close(pb)
    if (JoinAll == TRUE) {
        VCFs <- VCFs %>% bind_rows()
    }
    VCFs
}

#' @description
#' A function read a VCF file removing the header (denoted by the ##)
#'
#' @title process_folder
#' @param folder String containing the path to a folder of VCF files.
#' @return a dataframe of processed variants from all files in a folder.
#' @keywords VCF
#' @importFrom tibble remove_rownames
#' @export
#' @examples
#'
#' process_folder(file)

process_folder <- function(folder, pattern = ".annotated.vcf", ...) {
    VCFs <- list.files(folder, pattern = pattern, full.names = TRUE)
    process_VCFs(VCFs, ...)
}

#' A function to read a VCF file whilst removing the header (denoted by the ##)
#'
#' @title read_VCF
#' @param file The path to the file to be read
#' @param pass whether to filter out non PASS variants (default: TRUE)
#' @param pick whether to filter out the annotations that weren't picked (default: TRUE)
#' @param verbose print out extra text for debugging (verbose: FALSE)
#' @return A dataframe of the VCF skipping the commented lines
#' @keywords VCF
#' @importFrom readr read_delim
#' @importFrom dplyr filter rename
#' @export
#'
#' @examples
#'
#' read_VCF(file)

read_VCF <- function(file, pass = TRUE, verbose = FALSE, ...) {
    if (verbose) {
        print(paste0("Reading file ", file))
    }
    VCF <- read_delim(file,
                      delim = "\t",
                      ## This skips the comments at the top
                      comment = "##",
                      ## this stops it from printing out comments
                      show_col_types = FALSE) %>%
        rename("CHROM" = "#CHROM")
    if (pass) {
        if (verbose) {
            print("Filtering variants by PASS")
        }
        VCF <- VCF %>%
            filter(FILTER == "PASS")
    }
    VCF <- set_sample_name(VCF, file, ...)
    if (verbose) {
        print("Filtering large structural variants")
    }
    VCF <- VCF %>%
        filter(! grepl("<.*>", ALT))
    VCF
}

#' A function to process read counts from VCF files
#'
#' @title process_counts
#' @param file The path to the file to be read
#' @param VCF the VCF file with the tumour format split
#' @param Caller the caller the data is from. If empty determines from filepath (default: "")
#' @param verbose print out extra text for debugging (verbose: FALSE)
#' @return A dataframe of the VCF with ADR, RDR, DP and AF columns
#' @keywords VCF
#' @export
#'
#' @examples
#'
#' process_counts(file, VCF)

process_counts <- function(file, VCF, Caller = "", ...) {
    if (grepl("Strelka2", file) | Caller == "Strelka2") {
        VCF <- process_Strelka2(VCF, ...)
    } else if (grepl("Mutect2", file) | Caller == "Mutect2") {
        VCF <- process_Mutect2(VCF, ...)
    } else if (grepl("Varscan2", file) | Caller == "Varscan2") {
        VCF <- process_Varscan2(VCF, ...)
    } else if (grepl("VarDict", file) | Caller == "VarDict") {
        VCF <- process_VarDict(VCF, ...)
    } else {
        errorCondition("Cannot determine caller, stopping")
    }
    VCF
}

#' A function to split the tumour data from a VCF into columns
#'
#' @title split_format
#' @param VCF the VCF file as a dataframe
#' @return A dataframe of the tumour counts as columns
#' @keywords VCF
#' @importFrom dplyr filter bind_rows arrange
#' @importFrom stringr str_split_1
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#'
#' split_format(file, VCF)

split_format <- function(VCF, ...) {
    lapply(unique(VCF$FORMAT),function(format) {
        format_head <- str_split_1(format, ":")
        filter(VCF, FORMAT == format) %>%
            separate(TUMOR, into = format_head, ":", )}) %>%
        bind_rows() %>%
        arrange(CHROM, POS, REF, ALT)
}

#' A function to give an error for sample column finding
#'
#' @title split_format
#' @param colnames column names for error message
#' @param pattern pattern for error message
#' @return Error
#' @keywords VCF
#'

raise_sample_pattern_error <- function(colnames, pattern) {
    stop(cat("No column matches: ",
             pattern,
             "\nColumns: ",
             paste(colnames,
                   collapse = ", ")))
}

#' A function to determine the sample name and change to columns for the tumour and normal data to TUMOR and NORMAL if they aren't already
#'
#' @title set_sample_name
#' @param VCF the VCF file as a dataframe
#' @param filename the name of the VCF file
#' @param sample the sample ID. If NULL this will be determined by the filename (default: NULL)
#' @param tumourPattern The pattern for the tumour sample ids. Used to rename the tumour sample column to tumour (default: NULL)
#' @param normalPattern The pattern for the normal sample ids. Used to rename the normal sample column to NORMAL (default: NULL)
#' @param override_Varscan_handling Logical if TRUE will override additional varscan specific handelling
#' @param tumourOnly logical if TRUE no normal sample is present
#' @return A dataframe of the tumour counts as columns
#' @importFrom dplyr mutate rename
#' @keywords VCF
#' @export
#'
#' @examples
#'
#' set_sample_name(VCF, filename, tumourPattern = "tumour", normalPattern = "normal")
#'

set_sample_name <- function(VCF,
                            filename,
                            sample = NULL,
                            tumourPattern = NULL,
                            normalPattern = NULL,
                            override_Varscan_handling = FALSE,
                            tumourOnly = FALSE,
                            ...) {
    ##If there is no sample argument
    ##  take the file name and make it the Sample column
    if (is.null(sample)) {
        VCF <-
            dplyr::mutate(
                VCF,
                Sample = gsub("\\.[A-z]+[0-9]{0,1}\\.annotated.vcf.*", "", basename(filename)),
                .before = "INFO"
            )
        ##otherwise
        ##  make the sample argument the sample column
    } else {
        VCF <- dplyr::mutate(VCF, Sample = sample)
    }
    ## If TUMOR ISN'T a column name
    if (! "TUMOR" %in% colnames(VCF)) {
        ## IF Sample1 is in the VCF columns AND
        ## Sample1 ISN'T in the file name
        if ("Sample1" %in% colnames(VCF) &
            ! override_Varscan_handling &
            tumourOnly) {
            tumourCol <- "Sample1"
        }
        ## If there is a tumour pattern
        else if (! is.null(tumourPattern)) {
            ## if tumourPattern is in the colnames
            if (sum(grepl(tumourPattern, colnames(VCF))) > 0) {
                ## get the column that matches tumourPattern
                tumourCol <- colnames(VCF)[grepl(tumourPattern, colnames(VCF))]
                ## change the name of this column to TUMOR

            } else {
                raise_sample_pattern_error(colnames(VCF), tumourPattern)
            }
            ## if there isn't a tumour pattern
        } else {
            ## get the column containing the sample name
            tumourCol <- colnames(VCF)[grepl(unique(VCF$Sample), colnames(VCF))]
        }
        ## change the name of this column to TUMOR
        VCF <- rename(VCF, "TUMOR" = tumourCol)
    }
    ## if NORMAL isn't in colnames
    if (! "NORMAL" %in% colnames(VCF)) {
        ## if there is a normal pattern
        if (! is.null(normalPattern)) {
            if (sum(grepl(normalPattern, colnames(VCF)))>0) {
                ## get the column that matches the normalPattern
                normalCol <- colnames(VCF)[grepl(normalPattern, colnames(VCF))]
                ## change this column to NORMAL
                VCF <- rename(VCF, "NORMAL" = normalCol)
            }
            else {
                raise_sample_pattern_error(colnames(VCF), normalPattern)
            }
        }
    }
    VCF
}

#' A function to split VEP annotations into columns
#'
#' @title split_vep
#' @param VCF the VCF file as a dataframe
#' @param header Vector of the names for the new annotations columns. Must be the same length as the number of annotation columns.
#' @return A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
#' @keywords VCF
#' @importFrom dplyr rename bind_rows
#' @importFrom stringr str_split_1
#' @export
#'
#' @examples
#'
#' split_vep(VCF, header = FALSE)

split_vep <- function(VCF, header = FALSE, pick = TRUE, verbose = FALSE, ...) {
    if (header == FALSE) {
        ## This is the headings that should come out of the snakemake pipeline
        header = c("Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
                   "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
                   "cDNA_position","CDS_position", "Protein_position",
                   "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
                   "STRAND", "FLAGS", "PICK", "VARIANT_CLASS", "SYMBOL_SOURCE",
                   "HGNC_ID", "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL",
                   "TSL", "APPRIS", "CCDS", "ENSP","SWISSPROT", "TREMBL",
                   "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO", "SIFT", "PolyPhen",
                   "DOMAINS", "miRNA", "AF", "AFR_AF", "AMR_AF", "EAS_AF",
                   "EUR_AF", "SAS_AF", "gnomADe_AF", "gnomADe_AFR_AF",
                   "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
                   "gnomADe_FIN_AF", "gnomADe_NFE_AF", "gnomADe_OTH_AF",
                   "gnomADe_SAS_AF", "gnomADg_AF", "gnomADg_AFR_AF",
                   "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
                   "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF",
                   "gnomADg_NFE_AF", "gnomADg_OTH_AF", "gnomADg_SAS_AF",
                   "MAX_AF","MAX_AF_POPS", "CLIN_SIG", "SOMATIC", "PHENO",
                   "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS",
                   "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS", "LoFtool")
    }
    ## Iterate over the rows
    VCF <- apply(VCF, 1, function(Row) {
        ## Split INFO by |
        SplitRow <- Row[["INFO"]] %>% str_split_1(pattern = "\\|")
        ## Remove the columns you expect from the caller
        VEP <-
            SplitRow[-c(1:(length(SplitRow) %% length(header)))] %>%
            ## create a matrix with nrows the same length as the headers
            matrix(nrow = length(header)) %>%
            ## transpose this and make it a dataframe
            t() %>%
            as.data.frame() %>%
            ## set the names to those in the header
            setNames(header)
        ## bind the full row to the trimmed INFO column and VEP annotations
        cbind(Row %>% t(),
              INFO.Trim = paste(SplitRow[c(1:(length(SplitRow) %% length(header)))],
                                collapse = "|"),
              VEP)
    }) %>%
            bind_rows() %>%
        ## AF is the 1000 genomes total pop frequency, so rename this
        rename("AF_ALL" = "AF") %>%
        mutate(
            AF_ALL = ifelse(is.na(as.numeric(AF_ALL)), 0, as.numeric(AF_ALL)),
            gnomADe_AF = ifelse(is.na(as.numeric(gnomADe_AF)), 0, as.numeric(gnomADe_AF)),
            gnomADg_AF = ifelse(is.na(as.numeric(gnomADg_AF)), 0, as.numeric(gnomADg_AF))
        )
    if (pick) {
        if (verbose) {
            print("Filtering variants by PICK")
        }
        VCF <- VCF %>%
            filter(PICK == "1")
    }
    VCF
}

#' A function to split VEP annotations into columns
#'
#' @title process_AlleleMix
#' @param row A row of a VCF file
#' @param RA The allele to return (REF/ALT) Default = "ALT")
#' @return A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
#' @keywords VCF
#' @export
#'
#' @examples
#'
#' process_AlleleMix(row, RA = "ALT")

process_AlleleMix <- function(row, RA = "ALT", ...) {
    Allele = row[RA]
    Depth = row["DP"] %>% as.numeric()
    Counts = row[c("AU", "TU", "GU", "CU")] %>% gsub(",.*", "", .) %>%
        as.numeric() %>%
        setNames(c("A", "T", "G", "C"))
    Counts[Allele]
}

#' A function to split VEP annotations into columns
#'
#' @title process_Strelka2
#' @param VCF a dataframe of VCF records from Strelka with the TUMOR columns separated
#' @return A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
#' @keywords VCF
#' @importFrom dplyr mutate
#' @export
#'
#' @examples
#'
#' process_Strelka2(row)

process_Strelka2 <- function(VCF, ...) {
    ## Add allele frequencies
    VCF <- VCF %>% mutate(
        RDP = apply(VCF, 1, function(row){
            ifelse(is.na(row["TAR"]),
                   {process_AlleleMix(row, RA = "REF")},
                   row["TAR"] %>% gsub(",.*","", .) %>% as.numeric())
        }),
        ADP = apply(VCF, 1, function(row){
            ifelse(is.na(row["TIR"]),{
                process_AlleleMix(row, RA = "ALT")},
                row["TIR"] %>% gsub(",.*","", .) %>% as.numeric())
        }),
        DP = as.numeric(DP),
        AF = as.numeric(ADP)/as.numeric(DP)) %>%
        mutate(POS = as.numeric(POS))
    VCF
}

#' A function to split VEP annotations into columns
#'
#' @title process_Varscan2
#' @param VCF a dataframe of VCF records from Varscan2 with the TUMOR columns separated
#' @return A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
#' @keywords VCF
#' @importFrom dplyr mutate filter
#' @export
#'
#' @examples
#'
#' process_Varscan2(row)

process_Varscan2 <- function(VCF, sample = NULL, tumourOnly = FALSE, ...) {

    if (! tumourOnly ) {
        VCF <- VCF %>%
            ## SS=2 in INFO retains the SOMATIC variants
            filter(grepl("SS=2", INFO))
    }
    VCF <- VCF %>%
        mutate(RDP = as.numeric(RD),
               ADP = as.numeric(AD),
               DP = as.numeric(DP),
               AF = ADP/DP) %>%
        mutate(POS = as.numeric(POS))
    if (! tumourOnly) {
        VCF <- VCF %>%
            separate("DP4", c("RDF", "RDR", "ADF", "ADR"), convert = TRUE)
    } else if (tumourOnly) {
        VCF <- VCF %>%
            mutate(RDF = as.numeric(RDF),
                   RDR = as.numeric(RDR),
                   ADF = as.numeric(ADF),
                   ADR = as.numeric(ADR))
    }
    VCF
}

#' A function to split VEP annotations into columns
#'
#' @title process_Mutect2
#' @param VCF a dataframe of VCF records from Mutect2 with the TUMOR columns separated
#' @return A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
#' @keywords VCF
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#'
#' process_Mutect2(row)

process_Mutect2 <- function(VCF, ...) {
    VCF <- VCF %>%
        separate(AD, into = c("RDP", "ADP"), sep = ",") %>%
        mutate(DP = as.numeric(DP),
               RDP = as.numeric(RDP),
               ADP = as.numeric(ADP),
               AF = as.numeric(AF)) %>%
        mutate(POS = as.numeric(POS)) %>%
        separate("SB", c("RDF", "RDR", "ADF", "ADR"), convert = TRUE)
    VCF
}

#' A function to split VEP annotations into columns
#'
#' @title process_VarDict
#' @param VCF a dataframe of VCF records from VarDict with the TUMOR columns separated
#' @return A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
#' @keywords VCF
#' @importFrom dplyr mutate filter
#' @export
#'
#' @examples
#'
#' process_VarDict(row)

process_VarDict <- function(VCF, tumourOnly = FALSE, ...) {
    VCF <- VCF %>%
        separate(AD, into = c("RDP", "ADP"), sep = ",") %>%
        mutate(DP = as.numeric(DP),
               RDP = as.numeric(RDP),
               ADP = as.numeric(ADP),
               AF = as.numeric(AF)) %>%
        mutate(POS = as.numeric(POS)) %>%
        separate("RD", c("RDF", "RDR"), convert = TRUE) %>%
        separate("ALD", c("ADF", "ADR"), convert = TRUE)
    if ( ! tumourOnly ) {
        VCF <- filter(VCF, grepl("STATUS=.*Somatic", INFO))
    }
    VCF
}
