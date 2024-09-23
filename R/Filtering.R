#' A function to filter processed VCF calls
#'
#' @title filter_vars
#' @param VCF the VCF file process through the pipeline
#' @param biotype A vector of biotypes to select (default: c("protein_coding"))
#' @param impact A vector of biotypes to select (default: c("HIGH", "MODERATE"))
#' @param existing bool of whether to filter by dbSNP/COSMIC (default: FALSE)
#' @param population numeric value to filter 1000g/gnomAD (default: 0.01)
#' @param vaf numeric value to filter vaf (default: 0)
#' @param depth numeric value to filter depth (default: 1)
#' @param alt.depth numeric value to filter alt allele depth (default: 1)
#' @param returnDF bool return dataframe rather than vector (default: FALSE)
#' @param returnStats return a filtering summary (default: FALSE)
#' @param returnPlot return a plot summary of the filtering (default: FALSE)
#' @param plotStats print a plot summary of the filtering (default: FALSE)
#' @return A dataframe of the filtered VCF (or a plot/overview of the filtering)
#' @keywords VCF
#' @importFrom dplyr select arrange
#' @export
#' @examples
#'
#' filter_vars(VCF)
#' filter_vars(VCF, returnStats = FALSE)
#' filter_vars(VCF, af = 0.1)
#'

filter_vars <- function(VCF,
                        returnStats = FALSE,
                        returnPlot = FALSE,
                        plotStats = FALSE,
                        printStats = FALSE,
                        ...) {
    VCF <- add_filters(VCF, ...)
    Stats <- VCF %>%
        dplyr::select(starts_with("Filter."))
    if (printStats == TRUE) {
        print(summarise_filters(Stats, ...))
    }
    if (plotStats == TRUE) {
        print(plot_filters(Stats, ...))
    }
    if (returnStats == TRUE) {
        summarise_filters(Stats, ...)
    } else if (returnPlot == TRUE) {
        plot_filters(Stats, ...)
    } else {
        VCF[Stats %>% rowSums() == ncol(Stats), ] %>%
            dplyr::select(!starts_with("Filter.")) %>%
            arrange(CHROM, POS, REF, ALT)
    }
}

#' A function to return the chosen filters
#'
#' @title return_filters
#' @param biotype A vector of biotypes to select (default: c("protein_coding"))
#' @param impact A vector of biotypes to select (default: c("HIGH", "MODERATE"))
#' @param existing bool of wether to filter by dbSNP/COSMIC (default: FALSE)
#' @param population numeric value to filter 1000g/gnomAD (default: 0.01)
#' @param vaf numeric value to filter vaf (default: 0)
#' @param depth numeric value to filter depth (default: 1)
#' @param alt.depth numeric value to filter alt allele depth (default: 1)
#' @param returnDF bool return dataframe rather than vector (default: FALSE)
#' @return A list of the selected filters
#' @keywords VCF
#' @importFrom tibble rownames_to_column
#' @export
#' @examples
#'
#' return_filters()
#' return_filters(vaf = 0.10, depth = 50, alt.depth = 5)
#' return_filters(impact = c("LOW", "MODIFIER", "HIGH", "MODERATE"))
#'

return_filters <- function(biotype = c("protein_coding"),
                           impact = c("HIGH", "MODERATE"),
                           existing = TRUE,
                           population = 0.01,
                           vaf = 0,
                           depth = 1,
                           alt.depth = 1,
                           returnDF = FALSE,
                           ...) {
    filters <- list(
        "biotype" = biotype,
        "impact" = impact,
        "existing" = existing,
        "population" = population,
        "vaf" = vaf,
        "depth" = depth,
        "alt.depth" = alt.depth
    )
    if (returnDF == TRUE) {
        filters %>%
            lapply(paste, collapse = ", ") %>%
            unlist() %>%
            as.data.frame() %>%
            setNames("Value") %>%
            rownames_to_column("Filter")
    } else {
        filters
    }
}

#' A function to add columns for each of the filters with TRUE/FALSE if the variant passes/fails
#'
#' @title add_filters
#' @param VCF a dataframe of variants processed by the pipeline.
#' @return A dataframe of the VCF with TRUE/FALSE columns for each filter
#' @keywords VCF
#' @importFrom dplyr mutate
#' @export
#' @examples
#'
#' add_filters(VCF)
#'

add_filters <- function(VCF, ...) {
    filters = return_filters(...)
    VCF <- mutate(
        VCF,
        Filter.biotype = BIOTYPE %in% filters$biotype,
        Filter.impact = IMPACT %in% filters$impact,
        Filter.population = (
            (AF_ALL <= filters$population) &
                (gnomADe_AF <= filters$population) &
                (gnomADg_AF <= filters$population)
        ),
        Filter.vaf = AF >= filters$vaf,
        Filter.depth = DP >= filters$depth,
        Filter.alt.depth = ADP >= filters$alt.depth,
    )
    if (filters$existing == TRUE) {
        VCF <- mutate(VCF,
                      Filter.existing = (!grepl("rs", Existing_variation) |
                                             grepl("COS", ID)))
    }
    VCF
}

#' A function summarise the variants that pass the filters
#'
#' @title summarise_filters
#' @param Stats The boolean filter columns from the VCF
#' @return A list of the variants that pass each filter and how many variants pass how many filters
#' @keywords VCF
#' @importFrom dplyr mutate left_join
#' @importFrom tibble rownames_to_column
#' @export
#' @examples
#'
#' summarise_filters(VCF)
#'

summarise_filters <- function(Stats, ...) {
    list(
        Stats %>%
            colSums() %>%
            as.data.frame() %>%
            setNames("Passing") %>%
            rownames_to_column("Filter") %>%
            mutate(Filter = gsub("Filter.", "", Filter)) %>%
            left_join(return_filters(returnDF = TRUE, ...),
                      .,
                      by = "Filter"),
        table(Stats %>% rowSums) %>%
            as.data.frame() %>%
            setNames(c("FiltersPassed", "Variants"))
    ) %>% setNames(c("FilterSummary", "FiltersPassed"))
}

#' A function to plot the summary of the variants that pass the filters
#'
#' @title plot_filters
#' @param Stats The boolean filter columns from the VCF
#' @param title The title for the plot (default: "")
#' @return two barplots of the variants that pass each filter and the number of filters passed
#' @keywords VCF
#' @importFrom dplyr mutate left_join
#' @importFrom tibble rownames_to_column
#' @import patchwork
#' @export
#' @examples
#'
#' plot_filters(VCF)
#'

plot_filters <- function(Stats, title = "", ...) {
    PlotStats <- summarise_filters(Stats, ...)
    PlotStats$FilterSummary %>%
        ggplot(aes(x = Filter, y = Passing)) +
        geom_col() +
        geom_text(aes(
            x = Filter,
            y = Passing + (max(Passing) / 10),
            label = Passing
        )) +
        labs(title = title) +
        PlotStats$FiltersPassed %>%
        ggplot(aes(x = FiltersPassed, y = Variants)) +
        geom_col() +
        geom_text(aes(
            x = FiltersPassed,
            y = Variants + (max(Variants) / 10),
            label = Variants
        )) +
        plot_layout(ncol = 1)
}
