#' return_mean_pos_difference
#'
#' @title return_mean_pos_difference
#' @param df data frame of mutation positions
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange
#' @export
#'
#' @examples
#'
#' return_mean_pos_difference(df, idx)

return_mean_pos_difference <- function(df, ...) {
    positions <- df[, "POS"]
    differences <- positions[-1] - positions[-length(positions)]
    round(mean(differences))
}

#' return_cluster_format
#'
#' @title return_cluster_format
#' @param df data frame of mutation positions
#' @param idxs indexes to return
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange
#' @export
#'
#' @examples
#'
#' return_cluster_format(df, idx)


return_cluster_format <- function(df, idxs, ...){
    df <- df[idxs,]
    c(chr = unique(df$CHROM), start = min(df$POS), end = max(df$POS), num = nrow(df), mean_dist = return_mean_pos_difference(df))
}

#' mut_dist
#'
#' @title mut_dist
#' @param mutations data frame of mutation positions
#' @param distance max mean distance for variants (default: 1000)
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange
#' @export
#'
#' @examples
#'
#' mut_dist(mutations)


mut_dist <- function(mutations, distance = 1000, ...) {
    rv <- list()
    idx <- 1:6
    best.idx <- c()
    FINISHED = FALSE
    while (FINISHED == FALSE) {
        x <- mutations[idx, ]
        mean_diff <- return_mean_pos_difference(x)
        ## if the mean distance is below the threshold
        if (mean_diff <= distance) {
            best.idx <- idx
            ## check whether there are more mutations to add
            if (nrow(mutations) >= idx[length(idx)] + 1) {
                idx <- c(idx, idx[length(idx)] + 1)
            } else {
                rv[[length(rv) + 1]] <- return_cluster_format(mutations, best.idx)
                FINISHED = TRUE
            }
        } else {
            ## Need to store matches and start after matches
            if (!is.null(best.idx)) {
                rv[[length(rv) + 1]] <- return_cluster_format(mutations, best.idx)
                if (nrow(mutations) >= idx[length(idx)] + 5) {
                    idx <- (idx[length(idx)]):(idx[length(idx)] + 5)
                } else {
                    FINISHED = TRUE
                }
                ## no matches - increment by 1
            } else {
                if (nrow(mutations) >= idx[length(idx)] + 1) {
                    idx <- (idx[2]):(idx[length(idx)] + 1)
                } else {
                    FINISHED = TRUE
                }
            }
            best.idx <- c()
        }
    }
    bind_rows(rv)
}

#' process_mutation_distances
#'
#' @title process_mutation_distances
#' @param mutations data frame of mutation positions
#' @param obj object of class variant helpers
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange
#' @export
#'
#' @examples
#'
#' process_mutation_distances()
#'
#' process_mutation_distances(df)
#'
#' process_mutation_distances(list)
#'
#' process_mutation_distances(obj)

process_mutation_distances <- function(x) {
    UseMethod("process_mutation_distances")
}

#' @rdname process_mutation_distances
#' @export
process_mutation_distances.data.frame <- function(mutations, ...) {
    lapply(unique(mutations$Sample), function(sample) {
        lapply(unique(mutations$CHROM), function(chr) {
            muts <- mutations %>%
                dplyr::filter(Sample == sample) %>%
                dplyr::filter(CHROM == chr) %>%
                dplyr::select(CHROM, POS, REF, ALT) %>%
                dplyr::arrange(CHROM, POS) %>%
                dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)
            if(nrow(muts) > 6) {
                mut_dist(muts, ...) %>%
                    dplyr::mutate(Sample = sample)
            }
        })}) %>% bind_rows()
}

#' @rdname process_mutation_distances
#' @export
process_mutation_distances.list <- function(list, ...) {
    lapply(list, process_mutation_distances)
}

#' @rdname process_mutation_distances
#' @export
process_mutation_distances.VariantHelper <- function(obj, ...) {
    obj@misc$mutation_dist <- process_mutation_distances(obj@Callers)
    obj
}
