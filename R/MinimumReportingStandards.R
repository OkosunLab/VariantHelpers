#' Get minimum reporting stats from samtools
#'
#' @title get_samtools_stat
#' @param path path to QC folder
#' @param pattern pattern to find files (e.g. .con. to ignoire .raw.)
#' @param samtools_suffix folder with samtools output default: "/samtools_stats/"
#' @return data frame of metrics from samtools (metrics:
#' raw total sequences, reads mapped, reads unmapped, reads properly paired,
#' reads duplicated, average length, insert size average, pairs on different chromosomes"
#' cols:
#' Metrics, Value, Sample, Group)
#' @keywords VCF
#' @importFrom dplyr mutate filter select
#' @importFrom tidyr separate_wider_delim
#' @export
#'
#' @examples
#'
#' get_samtools_stat("QC", ".con.txt")

get_samtools_stat <- function(path, pattern, samtools_suffix = "/samtools_stats/", ...) {
    list.files(paste0(path, samtools_suffix), pattern = pattern, full.names = TRUE) %>%
        lapply(function(file) {
            readLines(file)[grep("^SN", readLines(file))] %>%
                as.data.frame() %>%
                setNames("Data") %>%
                tidyr::separate_wider_delim(Data,
                                     delim = "\t",
                                     names = c(NA,"Metric", "Value", NA),
                                     too_few = "align_start") %>%
                mutate(Metric = gsub(":", "", Metric)) %>%
                filter(Metric %in% c(
                    "insert size average", "average length", "pairs on different chromosomes",
                    "raw total sequences", "reads mapped", "reads properly paired", "reads unmapped",
                    "reads duplicated"
                ))%>%
                mutate(Value = as.numeric(Value),
                       Sample = gsub(paste0(".*\\/(.*)",pattern), "\\1", file))
        }) %>% bind_rows()
}

#' Get minimum reporting stats from picardtools
#'
#' @title get_picard_stat
#' @param path path to QC folder
#' @param picard_suffix folder with picardtools output default: /depth/picard/
#' @return data frame of metrics from picardtools (metrics:
#' MEAN_TARGET_COVERAGE, ZERO_CVG_TARGETS_PCT, PCT_TARGET_BASES_10X,
#' PCT_TARGET_BASES_30X, PCT_TARGET_BASES_100X;
#' cols:
#' Metrics, Value, Sample, Group)
#' @keywords VCF
#' @importFrom dplyr mutate
#' @export
#'
#' @examples
#'
#' get_picard_stat("QC")

get_picard_stat <- function(path, picard_suffix = "/depth/picard/", ...) {
    list.files(paste0(path, picard_suffix), full.names = TRUE) %>%
        lapply(function(file) {
            read_delim(file, delim = "\t", skip = 6, n_max = 1, show_col_types = FALSE) %>%
                dplyr::select(-BAIT_SET) %>%
                mutate(Sample = gsub(".*\\/(.*).txt", "\\1", file), .before = BAIT_TERRITORY) %>%
                pivot_longer(-1, names_to = "Metric", values_to = "Value") %>%
                filter(Metric %in% c("MEAN_TARGET_COVERAGE", "ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_100X"))
        }) %>%
        bind_rows()
}

#' Get the minimum reporting standards from picard and samtools
#'
#' @title get_min_reporting_stats
#' @param path path to QC folder
#' @param pattern pattern to find files (e.g. .con. to ignoire .raw.)
#' @param samtools_suffix folder with samtools output default: "/samtools_stats/"
#' @param picard_suffix folder with picard tools output default: /depth/picard/
#' @return data frame of metrics (cols: Metrics, Value, Sample, Group)
#' @keywords VCF
#' @importFrom dplyr mutate
#' @export
#'
#' @examples
#'
#' get_min_reporting_stats("QC")

get_min_reporting_stats <- function(...) {
    bind_rows(
        get_samtools_stat(...),
        get_picard_stat(...)
    ) %>%
        mutate(Value = as.numeric(Value),
               Metric =
                   case_when(
                       Metric == "MEAN_TARGET_COVERAGE" ~ "Mean target coverage",
                       Metric == "ZERO_CVG_TARGETS_PCT" ~ "Targets with 0 reads (%)",
                       Metric == "PCT_TARGET_BASES_10X" ~ "Target bases at 10X (%)",
                       Metric == "PCT_TARGET_BASES_30X" ~ "Target bases at 30X (%)",
                       Metric == "PCT_TARGET_BASES_100X" ~ "Target bases at 100X (%)",
                       .default = Metric
                   ) %>%
                   str_to_title() %>%
                   case_when(
                       . == "Insert Size Average" ~ "Insert Size Average (bp)",
                       . == "Average Length" ~ "Average Length (bp)",
                       .default = .
                   )
        ) %>%
        mutate(Group = case_when(
            grepl("Unmapped|Different|Duplicated", Metric) ~ "Failed Reads",
            grepl("^Reads|Sequences|Pairs", Metric) ~ "Reads",
            grepl("Length|Size", Metric) ~ "Fragment Size",
            grepl("Target Bases At", Metric) ~ "Depth",
            grepl("Targets with", Metric) ~ "Targets With 0 Reads",
            .default = Metric
        ) %>% factor(levels =
                         c("Reads", "Failed Reads", "Mean Target Coverage", "Depth",
                           "Fragment Size", "Targets With 0 Reads"))
        )
}


Plot_Colours = c(
    "Raw Total Sequences" = "olivedrab3",
    "Reads Mapped" = "olivedrab2",
    "Reads Properly Paired" = "olivedrab1",
    "Pairs On Different Chromosomes" = "orangered1",
    "Reads Duplicated" = "orangered3",
    "Reads Unmapped" = "orangered4",
    "Mean Target Coverage" = "coral",
    "Target Bases At 10x (%)" = "darkslategray4",
    "Target Bases At 30x (%)" = "darkslategray3",
    "Target Bases At 100x (%)" = "darkslategray2",
    "Insert Size Average (bp)" = "maroon3",
    "Average Length (bp)" = "orchid2",
    "Targets With 0 Reads (%)" = "brown3"
    )

#' Plot the minimum metrics
#'
#' @title plot_min_metrics
#' @param metrics dataframe of metrics produced by get_min_reporting_stats
#' @param tumourPattern optional pattern to split samples into tumour/normal default: NULL
#' @param overrideType optional ignore Type column whilst plotting
#' @param colours option vector of colours to use in the plot
#' @param independent string x, y or all. Vary scale within columns. Passed to facet_grid2
#' @return ggplot of supplied metrics
#' @keywords VCF
#' @importFrom dplyr mutate
#' @importFrom ggh4x facet_grid2
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' plot_min_metrics(metrics)
#' plot_min_metrics(metrics, tumourPattern = "FFPE")

plot_min_metrics <- function(metrics,
                             tumourPattern = NULL,
                             overrideType = FALSE,
                             colours = Plot_Colours,
                             ...) {
    if (! is.null(tumourPattern)) {
        metrics <- mutate(metrics,
                          Type = ifelse(grepl(tumourPattern, Sample),
                                        "Tumour",
                                        "Normal"))
    }
    plot <- ggplot(metrics,
                   aes(y = Sample,
                       x = Value,
                       fill = reorder(Metric,
                                      Value))) +
        geom_col(position = "dodge")
    if(! overrideType & "Type" %in% colnames(metrics)) {
        plot <- plot + ggh4x::facet_grid2(Type ~ Group, scales = "free", ...)
    } else {
        plot <- plot + ggh4x::facet_grid2( ~ Group, scales = "free_x", ...)
    }
    plot +
        theme_classic() +
        labs(title = "Minimum Reporting Metrics", fill = "", x = "") +
        scale_fill_manual(values = colours, breaks = names(colours))
}

