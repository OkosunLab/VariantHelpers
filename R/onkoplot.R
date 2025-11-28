#' @description
#' A collection of functions to plot the sample/gene level plots
#'
#' @rdname ancillary_plots
#' @title ancillary_plots
#' @param variants dataframe of variant calls
#' @param stat character matching the numerical column to plot (default: "VAF")
#' @param category column name to use for sample ID (default: Tumor_Sample_Barcode)
#' @import ggplot2
#' @examples
#'
#' plot_sample_stat(variants)
#'
plot_sample_stat <-
    function(variants,
             stat = "VAF",
             category = Tumor_Sample_Barcode,
             ...) {
        stat <- sym(stat)
        n_samp <- pull(variants, {{category}}) %>%
            levels() %>%
            length()
        plot <- variants %>%
            ggplot(aes(x = {
                {
                    category
                }
            }, y = !!stat)) +
            geom_vline(xintercept = seq(0.5, n_samp + 0.5), colour = "white") +
            geom_boxplot(outliers = FALSE) +
            geom_jitter() +
            scale_x_discrete(guide = "none", name = NULL, drop = FALSE) +
            theme(panel.grid.major.x = element_blank())
        plot
    }

#' @rdname ancillary_plots
#' @param variants dataframe of variant calls
#' @param Type column name to use for fill of plot (default: Variant_Classification)
#' @param COLOUR colours for the frequency plot (default: NULL)
#' @param category column name to use for sample ID (default: Tumor_Sample_Barcode)
#' @import ggplot2
#' @examples
#'
#' plot_sample_stat(variants)
#'
plot_sample_freq <-
    function(variants,
             category = Tumor_Sample_Barcode,
             Type = Variant_Classification,
             COLOUR = NULL,
             ...) {
        n_samp <- pull(variants, {{category}}) %>%
            levels() %>%
            length()
        plot <- variants %>%
            ggplot(aes(x = {
                {
                    category
                }
            }, fill = {
                {
                    Type
                }
            })) +
            geom_vline(xintercept = seq(0.5, n_samp + 0.5), colour = "white") +
            geom_bar() +
            scale_x_discrete(guide = "none", name = NULL, drop = FALSE) +
            theme(panel.grid.major.x = element_blank()) +
            labs(y = "TMB")
        if (is.null(COLOUR)) {
            plot <- plot +
                scale_fill_discrete(guide = "none", name = NULL,
                                    na.value = NA)
        } else {
            plot <- plot +
                scale_fill_manual(values = COLOUR, guide = "none",
                                  name = NULL, na.value = NA)
        }

        plot

    }

#' @rdname ancillary_plots
#' @param variants dataframe of variant calls
#' @param Type column name to use for fill of plot (default: Variant_Classification)
#' @param Gene column name to use for gene ID (default: Gene)
#' @import ggplot2
#' @examples
#'
#' plot_gene_freq(variants)
#'

plot_gene_freq <-
    function(variants,
             Gene = Gene,
             category = Tumor_Sample_Barcode,
             Type = Variant_Classification,
             COLOUR = NULL,
             ...) {
        Prevelence <- group_by(variants, {{Gene}}) %>%
            summarise(n = n()) %>%
            mutate(total = length(levels(pull(variants, {{category}}))),
                   percentage = round((n/total) * 100, 0) %>% paste0("%"))
        n_gene <- pull(variants, {{Gene}}) %>%
            unique() %>%
            length()
        plot <- variants %>%
            ggplot(aes(x = {{Gene}})) +
            geom_vline(xintercept = seq(0.5, n_gene + 0.5), colour = "white") +
            geom_bar(aes(fill = {{Type}})) +
            geom_text(data = Prevelence, aes(y = n + (total/2), label = percentage), hjust="inward") +
            coord_flip() +
            scale_x_discrete(guide = "none", name = NULL) +
            theme(panel.grid.major.y = element_blank()) +
            labs(y = "Frequency")
        if (is.null(COLOUR)) {
            plot <- plot +
                scale_fill_discrete(guide = "none", name = NULL)
        } else {
            plot <- plot +
                scale_fill_manual(values = COLOUR, guide = "none", name = NULL)
        }
        plot
    }

#' @rdname ancillary_plots
#' @param variants dataframe of variant calls
#' @param stat character matching the numerical column to plot (default: "VAF")
#' @param Gene column name to use for gene ID (default: Gene)
#' @import ggplot2
#' @examples
#'
#' plot_gene_stat(variants)
#'
plot_gene_stat <- function(variants,
                           stat = "VAF",
                           Gene = Gene,
                           ...) {
    stat <- sym(stat)
    n_gene <- pull(variants, {{Gene}}) %>%
        unique() %>%
        length()
    plot <- variants %>%
        ggplot(aes(x = !!stat,
                   y = {{Gene}})) +
        geom_hline(yintercept = seq(0.5, n_gene + 0.5), colour = "white") +
        geom_boxplot(outliers = FALSE) +
        geom_jitter() +
        scale_y_discrete(guide = "none", name = NULL) +
        theme(panel.grid.major.y = element_blank())
    plot
}

#' @rdname ancillary_plots
#' @import ggplot2
#' @param full.variants dataframe of unfiltered variants to get annotation information from
#' @param annotation_cols character/vector names of columns to add to the annotation plot (default: NULL)
#' @param annotation_colours named vector for annotation plot colours
#' @examples
#'
#' plot_gene_stat(variants, annotation_cols = c("sampleType", "donorID"))
#'
plot_annotation <- function(variants,
                            full.variants,
                            category = Tumor_Sample_Barcode,
                            annotation_cols = NULL,
                            annotation_colours = NULL,
                            ...) {
    if (! is.null(annotation_cols)) {
        n_samp <- pull(variants, {{category}}) %>%
            levels() %>%
            length()
        plot <- full.variants %>%
            mutate({{category}} := factor({{category}}, levels = levels(pull(variants, {{category}})))) %>%
            pivot_longer({{annotation_cols}}) %>%
            ggplot(aes(x = {{category}}, y = name, fill = value)) +
            geom_tile() +
            scale_x_discrete(guide = "none", name = "", drop = FALSE) +
            scale_y_discrete(name = "") +
            geom_vline(xintercept = seq(0.5, n_samp + 0.5), colour = "white") +
            theme(panel.grid.major.x = element_blank())
        if (! is.null(annotation_colours)) {
            plot <- plot +
                scale_fill_manual(values = annotation_colours, name = "")
        } else {
            plot <- plot +
                scale_fill_discrete(name = "")
        }
        plot <- list(plot)
    } else {
        plot <- list()
    }
    plot
}


#' @rdname multi_hit
#' @title set_multi_hit
#' @description
#' A function to convert multiple hits on the same gene in the same sample to a single entry labelled "multi hit"
#'
#' @param variants dataframe of variant calls
#' @param stat character matching the numerical column to plot (default: "VAF")
#' @param MultiHit bool collapse multi hits to multihit (default: TRUE)
#' @param MultiHitSG bool collapse multi hits to multihit in this plot. If NULL defaults to MultiHit (default: NULL)
#' @param Gene column name to use for gene ID (default: Gene)
#' @param Type column name to use for fill of plot (default: Variant_Classification)
#' @param category column name to use for sample ID (default: Tumor_Sample_Barcode)
#' @param multi_hit_label string to use for multi hits (default: multi_hit)
#' @import ggplot2
#' @examples
#'
#' set_multi_hit(variants)
#'

set_multi_hit <- function(variants,
                          category = Tumor_Sample_Barcode,
                          Gene = Gene,
                          Type = Variant_Classification,
                          multi_hit_label = "multi_hit",
                          ...) {
    variants <- variants %>%
        group_by({{Gene}}, {{category}}) %>%
        summarize(n = n(), {{Type}} := paste({{Type}}, sep = "|")) %>%
        mutate({{Type}} := ifelse(n == 1, {{Type}}, multi_hit_label)) %>%
        unique()
    variants
}

#' @title plot_gene_sample_mutations
#' @description
#' A function to generate a tile plot of the variants in each sample.
#'
#' @param variants dataframe of variant calls
#' @param stat character matching the numerical column to plot (default: "VAF")
#' @param MultiHit bool collapse multi hits to multihit (default: TRUE)
#' @param MultiHitSG bool collapse multi hits to multihit in this plot. If NULL defaults to MultiHit (default: NULL)
#' @param MultiHitS bool collapse multi hits to multihit in the sample plot. If NULL defaults to MultiHit (default: NULL)
#' @param MultiHitG bool collapse multi hits to multihit in the gene plot. If NULL defaults to MultiHit (default: NULL)
#' @param Gene column name to use for gene ID (default: Gene)
#' @param Type column name to use for fill of plot (default: Variant_Classification)
#' @param category column name to use for sample ID (default: Tumor_Sample_Barcode)
#' @import ggplot2
#' @examples
#'
#' plot_gene_stat(variants)
#'

plot_gene_sample_mutations <-
    function(variants,
             category = Tumor_Sample_Barcode,
             MultiHit = TRUE,
             MultiHitSG = NULL,
             Gene = Gene,
             Type = Variant_Classification,
             COLOUR = NULL,
             xlab = "Sample",
             xnames = FALSE,
             ...) {
        if ( is.null(MultiHitSG) ) {
            MultiHitSG = MultiHit
        }
        if ( MultiHitSG ) {
            variants <- set_multi_hit(variants,
                                      Type = {{Type}},
                                      category = {{category}},
                                      Gene = {{Gene}}, ...)
        }
        plot <- variants %>%
            mutate({{Type}} := factor({{Type}}, levels = names(COLOUR))) %>%
            ggplot(aes(x = {{category}},
                       y = {{Gene}},
                       fill = {{Type}}
                       )) +
            geom_tile(show.legend = TRUE)
        if (! is.null(COLOUR)) {
            plot <- plot +
                scale_fill_manual(values = COLOUR, breaks = names(COLOUR), drop = FALSE)
        }
        n_gene <- pull(variants, {{Gene}}) %>%
            unique() %>%
            length()
        n_samp <- pull(variants, {{category}}) %>%
            levels() %>%
            length()
        plot <- plot +
            theme(panel.grid = element_blank()) +
            geom_hline(yintercept = seq(0.5, n_gene + 0.5), colour = "white")  +
            geom_vline(xintercept = seq(0.5, n_samp + 0.5), colour = "white")
        if (! xnames) {
            plot <- plot +
                scale_x_discrete(guide = "none", drop = FALSE)
        } else {
            plot <- plot +
                scale_x_discrete(guide = guide_axis(angle = 90), drop = FALSE)
        }
        if ( ! is.null(xlab)) {
            plot <- plot +
                labs(x = xlab)
        }
        plot
    }



#' @title sample_gene_plots
#' @rdname sample_gene_plots
#' @description
#' A collection of functions lists of the top or side plots.
#'
#' @param variants dataframe of variant calls
#' @param full.variants dataframe of all variant calls
#' @param frequency bool whether to plot the frequency
#' @param sample_frequency bool whether to plot the sample frequency. Overwritten by frequency if NULL (default: NULL)
#' @param gene_frequency bool whether to plot the sample frequency. Overwritten by frequency if NULL (default: NULL)
#' @param stats vector/character of the stats to plot (default: c("VAF"))
#' @param sample_stats vector/character of the sample level stats to plot. Overwritten by stats if NULL (default: NULL)
#' @param gene_stats vector/character of the sample level stats to plot. Overwritten by stats if NULL (default: NULL)
#' @param filter.vars whether to filter to mutations in just the sample of interest.
#' @import ggplot2
#' @examples
#'
#' plot_gene_stat(variants)
#'

sample_plot <-
    function(variants,
             full.variants,
             frequency = TRUE,
             sample_frequency = NULL,
             stats = c("VAF"),
             sample_stats = NULL,
             MultiHitS = FALSE,
             MultiHit = TRUE,
             filter.vars = FALSE,
             ...) {
        if (! filter.vars) {
            variants <- set_sample_order_full(variants, full.variants, ...)
        }
        if ( is.null(MultiHitS) ) {
            MultiHitS = MultiHit
        }
        rv = list()
        if (is.null(sample_stats)) {
            sample_stats = stats
        }
        if (is.null(sample_frequency)) {
            sample_frequency = frequency
        }
        if (!is.null(sample_stats)) {
            rv <- c(rv, lapply(sample_stats, function(s) {
                plot_sample_stat(variants = variants, stat = s, ...)
            }))
        }
        if ( MultiHitS ) {
            variants <- set_multi_hit(variants,
                                      ...)
        }
        if (sample_frequency) {
            rv <- c(rv, list(plot_sample_freq(variants,
                                              ...)))
        }
        rv
    }

#' @rdname sample_gene_plots
#'
gene_plot <-
    function(variants,
             frequency = TRUE,
             gene_frequency = NULL,
             stats = c("VAF"),
             gene_stats = NULL,
             MultiHitG = NULL,
             MultiHit = TRUE,
             ...) {
        if ( is.null(MultiHitG) ) {
            MultiHitG = MultiHit
        }
        rv = list()
        if (is.null(gene_stats)) {
            gene_stats = stats
        }
        if (is.null(gene_frequency)) {
            gene_frequency = frequency
        }
        if ( MultiHitG ) {
            variants.freq <- set_multi_hit(variants, ...)
        } else {
            variants.freq <- variants
        }
        if (gene_frequency) {
            rv <- c(rv, list(plot_gene_freq(variants.freq, ...)))
        }
        if (!is.null(gene_stats)) {
            rv <- c(rv, lapply(gene_stats, function(s) {
                plot_gene_stat(variants = variants, stat = s, ...)
            }))
        }
        rv
    }


#' @title Onkoplot
#' @description
#' Create an onkoplot
#'
#' @param variants dataframe of variant calls
#' @param stat character matching the numerical column to plot (default: "VAF")
#' @param category column name to use for sample ID (default: Tumor_Sample_Barcode)
#' @param Gene column name to use for gene ID (default: Gene)
#' @param Type column name to use for fill of plot (default: Variant_Classification)
#' @param frequency bool whether to plot the frequency
#' @param sample_frequency bool whether to plot the sample frequency. Overwritten by frequency if NULL (default: NULL)
#' @param gene_frequency bool whether to plot the sample frequency. Overwritten by frequency if NULL (default: NULL)
#' @param stats vector/character of the stats to plot (default: c("VAF"))
#' @param sample_stats vector/character of the sample level stats to plot. Overwritten by stats if NULL (default: NULL)
#' @param gene_stats vector/character of the sample level stats to plot. Overwritten by stats if NULL (default: NULL)
#' @import ggplot2
#' @export
#' @examples
#'
#' plot_gene_stat(variants)
#'

onkoplot <- function(variants, ...) {
    full.variants <- variants
    logging(text = "ordering samples", ...)
    variants <- order_samples_and_genes(variants, ...)
    logging(text = "generating palette", ...)
    COLOUR <- set_colour_palette(full.variants, ...)
    logging(text = "generating sample plots", ...)
    sample <- sample_plot(variants, full.variants, COLOUR = COLOUR,  ...)
    logging(text = "generating gene plots", ...)
    gene <- gene_plot(variants, COLOUR = COLOUR, ...)
    logging(text = "generating oncoplots", ...)
    onco <- plot_gene_sample_mutations(variants, COLOUR = COLOUR, ...)
    logging(text = "generating annotation plots", ...)
    annot <- plot_annotation(variants, full.variants, ...)
    logging(text = "printing plot", ...)
    c(sample,
      list(onco),
      gene,
      annot) %>% wrap_plots() +
        plot_layout(guides = "collect",
                    design = calculate_design(sample, gene, annot, ...))
}


#' @title set_colour_palette
#' @description
#' define the colour palette
#'
#' @param colour vector of colours for the fill
#' @param variants dataframe of variant calls
#' @param Type column name to use for fill of plot (default: Variant_Classification)
#' @import ggplot2
#' @export
#' @examples
#'
#' (variants)
#'

set_colour_palette <- function(variants, colour = NULL, Type = Variant_Classification, multi_hit_label = "multi_hit", ...) {
    if (is.null(colour)) {
        rv <- c("black",
                scales::hue_pal(h = c(100,300))(pull(variants, {{Type}}) %>%
                                                    unique() %>%
                                                    length())) %>%
                    setNames(c(multi_hit_label,
                               pull(variants, {{Type}}) %>%unique()))
    } else {
        rv <- colour
    }
    rv
}

#' @title set_sample_order_full
#' @description
#' set the sample order for the full variant list
#'
#' @param variants dataframe of variant calls
#' @param full.variants dataframe of variant calls
#' @param category the column to use for the x axis
#' @import ggplot2
#' @export
#' @examples
#'
#' set_sample_order_full(variants, full.variants)
#'
set_sample_order_full <- function(variants, full.variants, category = Tumor_Sample_Barcode, ...) {
    full.variants %>%
        mutate(
            {{category}} := factor({{category}},
                                   pull(variants, {{category}}) %>% levels()
            )
        )
}

