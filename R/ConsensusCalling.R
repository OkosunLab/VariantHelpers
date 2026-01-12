#' A function to generate an upset plot of the consensus variants.
#'
#' @title upset_by_caller
#' @return upset plot of the variants called by the callers.
#' @param object object of class VariantHelper
#' @param set column in the calls for the set (default = Caller)
#' @param id column in the calls for the unique id (default = sam_var_id a newly constructed column made of varID and Sample)
#' @param category name of a column in the metadata (default: Sample)
#' @param colours a vector of colours to label the plot by with
#' @keywords VCF
#' @export
#' @examples
#'
#' upset_by_caller(obj)

upset_by_caller <- function(object,
                            set=Caller,
                            id = sam_var_id,
                            colour = Sample,
                            col_vec = NULL, ...) {
    combine_calls(object) %>%
        mutate(sam_var_id = paste(varID, Sample)) %>%
        filter(!is.na(Caller)) %>%
        plot_upset(set = {{set}},
                   id = {{id}},
                   colour = {{colour}},
                   col_vec = col_vec,
                   ...)
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
#' @param remove_cols character vector of columns to remove
#' @keywords VCF
#' @export
#' @examples
#'
#' add_consensus(obj)

add_consensus <- function(object, min_caller = 3, remove_cols = c(""), ...) {
    consensusCalls <- count_vars_by_caller_sample(object, ...) %>%
        filter(n >= min_caller)
    object@Consensus <-
        filter_by_consensus(object, consensus = consensusCalls) %>%
        dplyr::select(-any_of(
            c(
                .$FORMAT %>% unique() %>% lapply(str_split_1, ":") %>% unlist() %>% unique(),
                "QUAL","INFO", "INFO.Trim", "FORMAT", "NORMAL", "AF", "RDP", "ADP", "RDF", "RDR", "ADF","ADR",
                remove_cols
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
#' @title add_variant_summary
#' @return object of class VariantHelper
#' @param object object of class VariantHelper
#' @param category name of a column in the metadata (default: Sample)
#' @param stat the statistic to be used. Must exist in the dataframes for each caller (default: AF)
#' @param FUN function to use for the summary (default: mean)
#' @keywords VCF
#' @export
#' @examples
#'
#' add_variant_summary(obj)

add_variant_summary <- function(object, ...) {
    object@Consensus <- left_join(object@Consensus,
                                  variant_stat_summary(object, ...)
    )
    object
}

#' @rdname add_variant_summary
#' @export
add_summary <- function(object, ...) {
    warning("add_summary is depricated in favour of add_variant_summary and will be removed in a future release")
    object@Consensus <- left_join(object@Consensus,
                                  variant_stat_summary(object, ...)
    )
    object
}

#' A function to add variants called by n callers or more to the consensus slot of the object
#'
#' @title variant_stat_summary
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
#' variant_stat_summary(obj)

variant_stat_summary <- function(object, category = Sample, stat = AF, FUN = mean, ...){
    name = paste(deparse(substitute(FUN)), deparse(substitute(stat)), sep = "_")
    combine_calls(object) %>%
        group_by(varID, {{category}}) %>%
        summarise(stat = FUN({{stat}}, na.rm = TRUE),.groups = "keep") %>%
        dplyr::rename(!!name := "stat")
}

#' @rdname variant_stat_summary
#' @export
stat_summary <- function(object, category = Sample, stat = AF, FUN = mean, ...){
    warning("stat_summary is depricated in favour of variant_stat_summary and will be removed in a future release")
    name = paste(deparse(substitute(FUN)), deparse(substitute(stat)), sep = "_")
    combine_calls(object) %>%
        group_by(varID, {{category}}) %>%
        summarise(stat = FUN({{stat}}, na.rm = TRUE),.groups = "keep") %>%
        dplyr::rename(!!name := "stat")
}
