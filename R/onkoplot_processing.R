#' @title get_weights
#' @rdname sample_ordering
#' @description
#' weight numbers to allow cumulative ordering
#'
#' @param number the maximum number of items to weight
#' @import ggplot2
#' @examples
#'
#' get_weights(variants, 10)
#'
## weight numbers to allow cumulative ordering
## i.e. anything at first potision will always be weighted higher than the cumulative sum of all the values below
## sum(score) for i in group 1 but not in any groups > sum(score) in everything but group 1
get_weights <- function(number, ...) {lapply(1:number, function(i) {
    1/(2^i)
}) %>% unlist()}


#' @title get_n_top_genes
#' @rdname sample_ordering
#' @description
#' get the most mutated genes
#'
#' @param number the maximum number of genes to return
#' @param collapseCategory whether to collapse the counts at a sample level, i.e. number of genes mutated one time OR MORE per samples (default: TRUE)
#'
#' @import ggplot2
#' @examples
#'
#' get_weights(variants, 10)
#'
## get the top n genes with the weights needed for sorting
get_n_top_genes <- function(vars, number = 10, category = Tumor_Sample_Barcode, Gene = Gene, collapseCategory = TRUE, ...) {
    if (collapseCategory) {
        vars <- vars %>%
            dplyr::select({{category}}, {{Gene}}) %>%
            unique()
    }
    vars %>%
        group_by({{Gene}}) %>%
        summarise(n = n(), .groups = "keep") %>%
        arrange(-n) %>%
        head(n = number) %>%
        ungroup() %>%
        mutate(weights = get_weights(nrow(.), ...))
}

#' @title VariantProcessing
#' @rdname VariantProcessing
#' @description
#' weight numbers to allow cumulative ordering
#'
#' @param number the maximum number of items to weight
#' @import ggplot2
#' @examples
#'
#' get_weights(variants, 10)
#'

filter_variants_by_genes <- function(variants, genes, Gene = Gene, ...) {
    variants %>%
        filter({{Gene}} %in% pull(genes, {{Gene}}))
}
#' @title order_sample_by_weight
#' @rdname sample_ordering
#' @description
#' weight numbers to allow cumulative ordering
#'
#' @param variants dataframe of variants
#' @param genes dataframe of genes with gene name, number of occurances and weight
#' @param category column to use to stratify samples (default: Tumor_Sample_Barcode)
#' @param Gene column to use to stratify genes (default: Gene)
#' @import ggplot2
#' @examples
#'
#' get_weights(variants, 10)
#'
order_sample_by_weight <- function(variants, genes,
                                   category = Tumor_Sample_Barcode,
                                   Gene = Gene,
                                   reorder_by = NULL, ...) {
    variants.weight <- filter_variants_by_genes(variants, genes, Gene = {{Gene}}, ...) %>%
        dplyr::select({{category}}, {{Gene}}) %>%
        unique() %>%
        mutate(value = 1) %>%
        arrange(match({{Gene}}, pull(genes, {{Gene}}))) %>%
        pivot_wider(names_from = {{category}}, values_from = value, values_fill = 0) %>%
        column_to_rownames(colnames(.)[1])
    missingSamples <- unique(pull(variants, {{category}}))[
        !unique(pull(variants, {{category}})) %in% colnames(variants.weight)]
    if (length(missingSamples > 0)) {
        variants.weight <- bind_cols(variants.weight,
                                     lapply(missingSamples, function(sample) {
                                         rep(0, nrow(variants.weight)) %>%
                                             as.data.frame() %>%
                                             setNames(sample)
                                     }) %>%
                                         bind_cols())
    }
    if (! is.null(reorder_by)) {
        variants.weight <- rbind(variants.weight,
                                 Category = dplyr::select(variants, {{category}}, !!as.name(reorder_by)) %>%
                                     unique() %>%
                                     mutate(!! reorder_by := factor(!!as.name(reorder_by),
                                                                    levels = table(!!as.name(reorder_by)) %>%
                                                                        sort() %>%
                                                                        names()) %>%
                                                as.numeric()) %>%
                                     pull(var = !!as.name(reorder_by), name = {{category}}) %>%
                                     .[colnames(variants.weight)])
        genes <- rbind(genes,
                       list("Category", 0, 1))
    }
    (variants.weight * pull(genes, weights)) %>%
        colSums() %>%
        sort(decreasing = TRUE)
}

#' @title set_order
#' @rdname sample_ordering
#' @description
#' weight numbers to allow cumulative ordering
#'
#' @import ggplot2
#' @examples
#'
#' get_weights(variants, 10)
#'

set_order <- function(variants, genes, samples, category = Tumor_Sample_Barcode, Gene = Gene, reorder_by = NULL, ...) {
    variants <- variants %>%
        mutate({{category}} := factor({{category}}, levels = names(samples)),
               {{Gene}} := factor({{Gene}}, levels = rev(pull(genes, {{Gene}}))))
    variants
}

#' @title VariantProcessing
#' @rdname sample_ordering
#' @description
#' weight numbers to allow cumulative ordering
#'
#' @import ggplot2
#' @examples
#'
#' get_weights(variants, 10)
#'

order_samples_and_genes <- function(variants, ...) {
    genes <- get_n_top_genes(variants, ...)
    samples <- order_sample_by_weight(variants, genes, ...)
    variants <- set_order(variants, genes, samples, ...)
    variants <- filter_variants_by_genes(variants, genes, ...)

}

#' @title calculate_design
#' @description
#' Work out the design string from the given top and side plots
#'
#' @param top the list of plots for the top of the plot
#' @param side the list of plots for the side of the plot
#' @param main.scale the relative size of the main plot (default: 6)
#' @import ggplot2
#' @examples
#'
#' calculate_design(top, side, , 10)
#'

calculate_design <- function(top, side, annot, main.scale = 6, ...) {
    c(
        lapply(LETTERS[1:length(top)], rep, main.scale) %>% lapply(paste, collapse = "") %>%
            paste(paste(rep(
                "#", length(side)
            ), collapse = ""), sep = ""),
        lapply(1:main.scale, function(i) {
            c(rep(LETTERS[length(top) + 1], main.scale),
              LETTERS[(length(top) + 2):(length(top) + 1 + length(side))]) %>%
                paste(collapse = "")
        }) %>% unlist(),
        rep(rep(LETTERS[length(top) + 2 + length(side)], main.scale), length(annot)) %>%
            paste(collapse = "") %>%
            paste(rep(paste(rep("#", length(side)), collapse = ""), length(annot)), sep = "")
    ) %>%
        paste(collapse = "\n")
}


#' @title logging
#'
#' @param verbose bool should the message be printed
#' @param text the text to plot
#'
#'
logging <- function(verbose = FALSE, text = "", ...) {
    if (verbose) {
        message(text)
    }
}
