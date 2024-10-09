#' A function to load the depth information from the Okosun lab alignment snakemake pipelines (or anything that uses GATK DepthOfCoverage)
#'
#' @title filter_vars
#' @param folders A vector/string of folders to look in for files (default: QC/depth/).
#' @param meta.data a dataframe containing aditional sample level information.
#' @param pattern a regex for finding the right files (default: .*sample_summary)
#' @return A dataframe of the GATK depth summaries.
#' @keywords VCF
#' @importFrom dplyr bind_rows filter left_join
#' @export
#' @examples
#'
#' load_depth_qc()
#'
#' load_depth_qc(folders = "WES/QC/depth", meta.data = MData)
#'

load_depth_qc <-
    function(folders = "QC/depth/",
             pattern = ".*sample_summary",
             meta.data = NULL,
             ...) {
        Files <-
            list.files(path = folders,
                       pattern = pattern,
                       full.names = TRUE)
        depth <- lapply(Files, read.csv) %>%
            bind_rows() %>%
            filter(sample_id != "Total")
        if (!is.null(meta.data)) {
            depth <- left_join(depth, meta.data)
        }
        depth
    }

#' A function to load the depth information from the Okosun lab alignment snakemake pipelines (or anything that uses GATK DepthOfCoverage)
#'
#' @title filter_vars
#' @param depth dataframe of depth of coverage statistics. If NULL will try to load data using load_depth_qc (default: NULL)
#' @param category column in depth to use for plotting
#' @return A ggplot object
#' @keywords VCF
#' @import ggplot
#' @import patchwork
#' @export
#' @examples
#'
#' plot_depth_qc(depth)
#'

plot_depth_qc <- function(depth = NULL, category = Sample, ...) {
    if (is.null(depth)) {
        load_depth_qc(...)
    }
    ## Density
    ggplot(depth,aes(x = mean)) +
        geom_density(alpha = 0.5) +
        labs(x = "Mean Base Coverage", fill = "Origin")
    ggplot(Depth,aes(y = mean, x = Type, fill = Type)) +
        geom_boxplot(outliers = FALSE) +
        geom_jitter() +
        labs(y = "Mean Base Coverage", fill = "Origin")
}

 +
    plot_layout(guides = "collect", ncol = 1) +
    expand_limits(y = 0) +
    ggplot(Depth,aes(y = mean, x = Patient, fill = Patient)) +
    geom_col() +
    facet_wrap(~factor(Type, levels=c("CSF", "PB", "Tumour", "Normal")),
               scales = "free_y", ncol = 1) +
    scale_fill_manual(values = ClinicalColours$Patient) +
    plot_layout(heights = c(1,1,4)) +
    scale_x_discrete(guide = guide_axis(angle = 45))
