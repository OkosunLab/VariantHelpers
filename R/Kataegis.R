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
    process_mutation_distances(obj@Callers)
}

#' mutation_dist_to_grange
#'
#' @title mutation_dist_to_grange
#' @param df of mutations distances
#' @param list of mutation distances dataframes
#' @param obj of type VariantHelper
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange rename
#' @export
#'
#' @examples
#'
#' mutation_dist_to_grange()
#'
#' mutation_dist_to_grange(df)
#'
#' mutation_dist_to_grange(list)
#'
#' mutation_dist_to_grange(obj)

mutation_dist_to_grange <- function(x) {
    UseMethod("mutation_dist_to_grange")
}
#' @rdname mutation_dist_to_grange
#' @export
mutation_dist_to_grange.data.frame <- function(df) {
    df %>%
        GenomicRanges::makeGRangesFromDataFrame(
            keep.extra.columns = TRUE,
            ignore.strand = TRUE,
            start.field = "start",
            end.field = "end")
}
#' @rdname mutation_dist_to_grange
#' @export
mutation_dist_to_grange.list <- function(list) {
    lapply(list, mutation_dist_to_grange)
}
#' @rdname mutation_dist_to_grange
#' @export
mutation_dist_to_grange.VariantHelper <- function(obj) {
    lapply(obj@misc$mutation_dist, mutation_dist_to_grange)
}


#' annotate_mutation_dist
#'
#' @title annotate_mutation_dist
#' @param mut_dist data frame of mutation positions
#' @param Genes GRanges object of genes to annotate
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange rename
#' @export
#'
#' @examples
#'
#' annotate_mutation_dist()
#'
#' annotate_mutation_dist(df)
#'
#' annotate_mutation_dist(list)
#'
#' annotate_mutation_dist(obj)


annotate_mutation_dist <- function(x) {
    UseMethod("annotate_mutation_dist")
}

#' @rdname annotate_mutation_dist
#' @export
annotate_mutation_dist.data.frame <- function(mut_dist, Genes = NULL, ...) {
    mut_dist <- mutation_dist_to_grange(mut_dist)
    if ( is.null(Genes) ) {
        Genes <- get_all_ens_genes(...)
    }
    Match <- GenomicRanges::findOverlaps(mut_dist, Genes) %>% as.data.frame()
    cbind(mut_dist[Match$queryHits,] %>% as.data.frame(),
          Genes[Match$subjectHits,] %>% as.data.frame() %>% dplyr::select(-c(1:5))) %>%
        unite("Locus", c(seqnames, band), sep = "") %>%
        group_by(Sample, start, end, num, mean_dist) %>%
        summarize(Locus = Locus %>% unique() %>% stringi::stri_remove_empty() %>% paste(collapse = ", "),
                  Gene = Gene %>% unique() %>% stringi::stri_remove_empty() %>% paste(collapse = ", "),
                  .groups = "drop") %>%
        dplyr::select(all_of(c("Sample", "Locus", "start", "end", "num", "mean_dist", "Gene"))) %>%
        dplyr::rename("Start Position (bp)" = start,
                      "End Position (bp)" = end,
                      "Number of Mutations" = num,
                      "Mean Intermutation Distance (bp)" = mean_dist)
}

#' @rdname annotate_mutation_dist
#' @export
annotate_mutation_dist.list <- function(list, ...) {
    lapply(list, annotate_mutation_dist, ...)
}

#' @rdname annotate_mutation_dist
#' @export
annotate_mutation_dist.VariantHelper <- function(obj, ...) {
    lapply(obj@misc$mutation_dist, annotate_mutation_dist, ...)
}

#' detect_kataegis
#'
#' @title detect_kataegis
#' @param df of mutations
#' @param list of mutation dataframes
#' @param obj of type VariantHelper
#' @param Genes GRanges object of genes to annotate
#' @keywords VCF
#' @importFrom dplyr mutate filter select arrange
#' @export
#'
#' @examples
#'
#' detect_kataegis(df)
#'
#' detect_kataegis(list)
#'
#' detect_kataegis(obj)

detect_kataegis <- function(x){
    UseMethod("detect_kataegis")
}

#' @rdname detect_kataegis
#' @export
detect_kataegis.data.frame <- function(df) {
    process_mutation_distances(df) %>%
        annotate_mutation_dist()
}

#' @rdname detect_kataegis
#' @export
detect_kataegis.list <- function(list) {
    process_mutation_distances(list) %>%
        annotate_mutation_dist()
}

#' @rdname detect_kataegis
#' @export
detect_kataegis.VariantHelper <- function(obj) {
    process_mutation_distances(obj) %>%
        annotate_mutation_dist()
}


