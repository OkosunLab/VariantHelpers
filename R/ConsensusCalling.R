#' A function to generate an upset plot of the consensus variants.
#'
#' @title upset_by_caller
#' @return upset plot of the variants called by the callers.
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @param colours a vector of colours to label the plot by with
#' @keywords VCF
#' @export
#' @examples
#'
#' upset_by_caller(obj)

upset_by_caller <- function(object,
                            category = Sample,
                            colours = NULL, ...) {
    if (is.null(colours)) {
        s.colour <- scale_fill_discrete(guide = "none")
    } else {
        s.colour <- scale_fill_manual(values = colours, guide = "none")
    }
    combine_calls(object, ...) %>%
        dplyr::select(varID, Caller, Sample, {{category}}, Consequence) %>%
        mutate(values = as.integer(1)) %>%
        pivot_wider(
            names_from = Caller,
            values_from = values,
            values_fill = as.integer(0)
        ) %>%
        as.data.frame() %>%
        ComplexUpset::upset(
            names(object@Callers),
            name = "caller",
            width_ratio = 0.2,
            base_annotations = list(
                'Intersection size' = intersection_size(
                    mapping = aes(fill = {{category}})
                ) +
                    s.colour
            ),
            set_sizes = (
                upset_set_size(geom = geom_bar(aes(
                    fill = {{category}}, x = group
                ),
                width = 0.8),) +
                    s.colour
            ))
}

#' A function to count the number of callers that called each variant.
#'
#' @title count_vars_by_caller_sample
#' @return dataframe of the number of times a variant was called
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @keywords VCF
#' @export
#' @examples
#'
#' count_vars_by_caller_sample(obj)

count_vars_by_caller_sample <- function(object, category = Sample, ...) {
    combine_calls(object, ...) %>%
        group_by(varID, {{category}}) %>%
        summarise(n = n(), .groups = "keep") %>%
        mutate(matchID = paste(varID, {{category}}))
}


#' A function to add variants called by n callers or more to the consensus slot of the object
#'
#' @title add_consensus
#' @return object of class VariantHelper with the consensus calls in the consensus slot.
#' @param object object of class VariantHelper
#' @param min_caller minimum number of callers needed to retain variant (default: 3)
#' @keywords VCF
#' @export
#' @examples
#'
#' add_consensus(obj)

add_consensus <- function(object, min_caller = 3, ...) {
    consensusCalls <- count_vars_by_caller_sample(object, ...) %>%
        filter(n >= min_caller)
    object@Consensus <-
        filter_by_consensus(object, consensus = consensusCalls) %>%
        dplyr::select(-any_of(
            c(
                .$FORMAT %>% unique() %>% lapply(str_split_1, ":") %>% unlist() %>% unique(),
                "QUAL","INFO", "INFO.Trim", "FORMAT", "NORMAL", "AF", "RDP", "ADP"
            )
        )) %>%
        mutate(value = 1) %>%
        pivot_wider(
            names_from = Caller,
            values_from = value,
            values_fill = 0
        )
    object@CallIDs <-
        paste(object@Consensus$varID, object@Consensus$Sample, sep = "_")
    object@Callers <- lapply(object@Callers, function(caller, category = Sample) {
        caller %>%
            unite("varID", c(CHROM, POS, REF, ALT), remove = FALSE) %>%
            filter(paste(varID, {{category}}, sep = "_") %in% object@CallIDs)
    })
    object
}

#' A function to retain only the varians in consensus
#'
#' @title filter_by_consensus
#' @return dataframe of the variants that exist in consensus
#' @param object object of class VariantHelper
#' @param consensus vector of Variant IDs for consensus variants
#' @keywords VCF
#' @export
#' @examples
#'
#' filter_by_consensus(obj)

filter_by_consensus <- function(object, consensus, category = Sample, ...) {
    combine_calls(object, ...) %>%
        filter(paste(varID, {{category}}) %in% consensus$matchID)
}

#' A function return the consensus calls from the object with the meta.data added.
#'
#' @title return_consensus
#' @return dataframe of consensus calls with the data in meta.data attached.
#' @param object object of class VariantHelper
#' @keywords VCF
#' @export
#' @examples
#'
#' return_consensus(obj)

return_consensus <- function(object) {
    left_join(object@meta.data, object@Consensus)
}

#' A function to add a summary of a statistic from the individual callers to the consensus dataframe
#'
#' @title add_summary
#' @return object of class VariantHelper
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @param stat the statistic to be used. Must exist in the dataframes for each caller (default: AF)
#' @param FUN function to use for the summary (default: mean)
#' @keywords VCF
#' @export
#' @examples
#'
#' add_summary(obj)


add_summary <- function(object, ...) {
    object@Consensus <- left_join(object@Consensus,
                                  stat_summary(object, ...)
    )
    object
}

#' A function to add variants called by n callers or more to the consensus slot of the object
#'
#' @title stat_summary
#' @return A dataframe of the consensus variant calls with the stat summary attached.
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @param stat the statistic to be used. Must exist in the dataframes for each caller (default: AF)
#' @param FUN function to use for the summary (default: mean)
#' @param min_caller minimum number of callers needed to retain variant (default: 3)
#' @keywords VCF
#' @importFrom dplyr rename
#' @export
#' @examples
#'
#' stat_summary(obj)

stat_summary <- function(object, category = Sample, stat = AF, FUN = mean, ...){
    name = paste(deparse(substitute(FUN)), deparse(substitute(stat)), sep = "_")
    combine_calls(object) %>%
        group_by(varID, {{category}}) %>%
        summarise(stat = FUN({{stat}}),.groups = "keep") %>%
        dplyr::rename(!!name := "stat")
}
