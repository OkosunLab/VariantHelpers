#' Hold variant call data
#'
#' This class will hold your vcf data
#' @title VariantHelper Class
#' @keywords VCF variants
#' @export
#' @examples
#' Community()

VariantHelper <- setClass(
    "VariantHelper",
    slots = c(
        Callers = "list",
        Samples = "vector",
        Consensus = "data.frame",
        CallIDs = "vector",
        meta.data = "data.frame",
        misc = "list",
        info = "list")
)

setMethod("show", "VariantHelper",
          function(object) {
              cat("Object of class VariantHelper\n",
                  "Contains calls from the following variant callers:\n",
                  names(object@Callers),
                  "\n",
                  "Samples: ", length(object@Samples))
          })


#' A function to generate an object of class VariantHelper for a vector of paths to folders
#'
#' @title variant_helper_from_folder_paths
#' @param paths a vector of paths to folders containing VCF files
#' @return object of class VariantHelper with the VCF files loaded
#' @keywords VCF
#' @export
#' @examples
#'
#' variant_helper_from_folder_paths(Paths)

variant_helper_from_folder_paths <- function(paths, ...) {
    Calls <- lapply(paths, function(path) {
        process_folder(path)
    }) %>%
        setNames(names(paths))
    rv <- variant_helper_from_list_of_calls(Calls, ...)
    return(rv)
}

#' A function to generate an object of class VariantHelper for a vector of paths to folders
#'
#' @title variant_helper_from_list_of_calls
#' @param Calls a list of processed VCF files
#' @param meta.data dataframe of metadata for samples (default: NULL)
#' @return object of class VariantHelper with the VCF info from Calls
#' @keywords VCF
#' @export
#' @examples
#'
#' variant_helper_from_list_of_calls(Calls)

variant_helper_from_list_of_calls <- function(Calls, meta.data = NULL, ...) {
    rv <- VariantHelper()
    rv@Callers <- Calls
    rv <- update_stats(rv)
    if (! is.null(meta.data)) {
        rv <- add_metadata(rv, metadata = meta.data)
    }
    return(rv)
}

#' A function to update the stats held by the object
#'
#' @title update_stats
#' @return object of class VariantHelper with the stats updated
#' @keywords VCF
#' @export
#' @examples
#'
#' update_stats(obj)

update_stats <- function(object, ...) {
    object@Samples <- lapply(object@Callers, function(caller){
        unique(caller$Sample)
    }) %>% unlist() %>% unique()
    return(object)
}

#' A wrapper to filter all VCFs loaded in obj@Callers
#'
#' @title update_stats
#' @param biotype A vector of biotypes to select (default: c("protein_coding"))
#' @param impact A vector of biotypes to select (default: c("HIGH", "MODERATE"))
#' @param existing bool of whether to filter by dbSNP/COSMIC (default: TRUE)
#' @param population numeric value to filter 1000g/gnomAD (default: 0.01)
#' @param vaf numeric value to filter vaf (default: 0)
#' @param depth numeric value to filter depth (default: 1)
#' @param alt.depth numeric value to filter alt allele depth (default: 1)
#' @return object of class VariantHelper with filters applied
#' @keywords VCF
#' @export
#' @examples
#'
#' filter_variants(obj, vaf = 0.5)

filter_variants <- function(object, ...) {
    object@Callers <- lapply(object@Callers, function(obj) {
        filter_vars(obj, ...)
    })
    object <- update_stats(object, ...)
    return(object)
}

#' A function to update the stats held by the object
#'
#' @title add_metadata
#' @return object of class VariantHelper with metadata added
#' @param metadata dataframe of sample metadata
#' @keywords VCF
#' @export
#' @examples
#'
#' add_metadata(obj, metadata)

add_metadata <- function(object, metadata, ...) {
    if (sum(object@Samples %in% metadata$Sample) == length(object@Samples)) {
        object@meta.data <- metadata
    } else {
        errorCondition("Samples in metadata don't match samples in object")
    }
    object
}

#' A function to return the data frames holding variant calls.
#' Can optionally return a subset of the callers defined by a vector
#'
#' @title return_callers
#' @return List of dataframes containing variant call information from all/selected callers.
#' @param object object of class VariantHelper
#' @param callers vector of the names of variant callers stored in the Callers slot.
#' @keywords VCF
#' @export
#' @examples
#'
#' return_callers(obj, metadata)

return_callers <- function(object, callers = NULL, ...) {
    if (is.null(callers)) {
        return(object@Callers)
    } else {
        return(object@Callers[callers])
    }
}

#' A function to take data from selected callers and combine them into a long dataframe
#'
#' @title combine_calls
#' @return Long data frame of calls from all the callers.
#' @param object object of class VariantHelper
#' @keywords VCF
#' @export
#' @examples
#'
#' combine_calls(obj)
#'
#' combine_calls(obj)

combine_calls <- function(object, ...) {
    return_callers(object, ...) %>%
        bind_rows(.id = "Caller") %>%
        left_join(object@meta.data, ., by = "Sample") %>%
        unite("varID", c(CHROM, POS, REF, ALT), remove = FALSE)
}

#' A function to plot the number of variants caller per sample (or other metadata column) per caller.
#'
#' @title variants_by_caller
#' @return ggplot of the number of variants.
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @keywords VCF
#' @export
#' @examples
#'
#' variants_by_caller(obj)

variants_by_caller <- function(object, category = Sample, ...) {
    combine_calls(object, ...) %>%
        ggplot(aes(x = {{category}}, fill = {{category}})) +
        geom_bar() +
        facet_wrap(~ Caller) +
        scale_x_discrete(guide = guide_axis(angle = 90))
}

#' A function to plot the chosen statistic for each variant per sample (or other metadata column) per caller.
#'
#' @title stats_by_caller
#' @return ggplot of the chose statistic of variants.
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @keywords VCF
#' @export
#' @examples
#'
#' stats_by_caller(obj)

stats_by_caller <- function(object, Stat = AF,  category = Sample, ...) {
    combine_calls(object, ...) %>%
        ggplot(aes(x = {{Stat}}, y = {{category}},
                   colour = {{category}}, fill = {{category}})) +
        ggridges::geom_density_ridges(alpha = 0.25) +
        facet_wrap(~Caller) +
        theme_classic()
}
