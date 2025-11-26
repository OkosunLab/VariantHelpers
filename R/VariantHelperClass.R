#' Hold variant call data
#'
#' This class will hold your vcf data
#' @title VariantHelper Class
#' @keywords VCF variants
#' @slot Callers list of dataframes containing call data
#' @slot Samples vecotr of sample names
#' @slot Consensus dataframe of consensus calls
#' @slot CallIDs vector
#' @slot meta.data dataframe of metadata
#' @slot misc list
#' @slot info list
#' @export
#' @examples
#'
#' VariantHelper()

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
        process_folder(path, ...)
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
    } else {
        rv <- add_metadata(rv,
                           metadata = as.data.frame(list(Sample = rv@Samples)))
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

#' A function to add metadata to the object
#'
#' @title add_variant_metadata
#' @return object of class VariantHelper with metadata added
#' @param metadata dataframe of sample metadata
#' @keywords VCF
#' @export
#' @examples
#'
#' add_variant_metadata(obj, metadata)

add_variant_metadata <- function(object, metadata, ...) {
    if (sum(object@Samples %in% metadata$Sample) == length(object@Samples)) {
        object@meta.data <- metadata
    } else {
        errorCondition("Samples in metadata don't match samples in object")
    }
    object
}

## Remove eventually to stop conflicts
#' @rdname add_variant_metadata
#' @export

add_metadata <- function(object, metadata, ...) {
    warning("add_metadata is depricated in favour of add_variant_metadata and will be removed in a future release")
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
#' @import ggplot2
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
#' @import ggplot2
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

#' A function to combine two instances of VariantHelper
#'
#' @title combine_objects
#' @return combined object of class VariantHelper
#' @param obj object of class VariantHelper
#' @param obj2 object of class VariantHelper
#' @keywords VCF
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#'
#' combine_objects(obj, obj2)

combine_objects <- function(obj, obj2) {
    rv <- obj
    rv@Consensus <- bind_rows(obj@Consensus,
                              obj2@Consensus)
    if (sum(duplicated(rv@Consensus)) > 0) {
        stop("Duplicated variants in Consensus")
    }
    callers <- unique(c(names(obj@Callers), names(obj2@Callers)))
    rv@Callers <- lapply(callers, function(caller) {
        # if (caller %in% names(obj@Callers)) {
        #
        # }
        rv <- dplyr::bind_rows(obj@Callers[[caller]],
                        obj2@Callers[[caller]])
        if (sum(duplicated(rv)) > 0) {
            stop(paste("Duplicated variants in", caller))
        }
        rv
    }) %>% setNames(callers)
    rv@meta.data <- dplyr::bind_rows(obj@meta.data,
                              obj2@meta.data)
    if (sum(duplicated(rv@meta.data)) > 0) {
        stop("Duplicated variants in Consensus")
    }
    rv <- update_stats(rv)
    rv
}
