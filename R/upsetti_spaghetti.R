#' upsetti_spaghetti
#'
#' This class will hold upset data
#' @name upsetti_spaghetti
#' @title upsetti_spaghetti
#' @keywords VCF variants upset
#' @slot plots a list of plot objects
#' @slot data a data.frame of the long form data
#' @export
#' @examples
#'
#' upsetti_spaghetti()

upsetti_spaghetti <- setClass(
    "upsetti_spaghetti",
    slots = c(
        plots = "list",
        data = "data.frame"
    )
)

setMethod("show", "upsetti_spaghetti",
          function(object) {
              if ("final" %in% names(object@plots)) {
                  print(object@plots$final)
              } else {
                  print("This spaghetti isn't very upsetti")
              }
          })

#' @rdname upsetti_spaghetti
#' @importFrom dplyr filter
#' @export

filter.upsetti_spaghetti <- function(object, ...) {
    dplyr::filter(object@data, ...)
}

#' @rdname upsetti_spaghetti
#' @importFrom dplyr mutate
#' @export
mutate.upsetti_spaghetti <- function(object, ...) {
    dplyr::mutate(object@data, ...)
}


#' A function to generate and process data into an upsetti_spaghetti class
#'
#' @title generate_upsetti
#' @return a new upsetti_spaghetti object
#' @param df a wide data.frame of data to plot
#' @param cols a vector of columns to pivot on
#' @keywords VCF
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter
#' @export
#' @examples
#'
#' generate_upsetti(df, cols = c("col1", "col2"))

generate_upsetti <- function(df, cols, ...) {
    rv <- upsetti_spaghetti()
    rv@data <- df %>%
        tidyr::pivot_longer(all_of(cols)) %>%
        dplyr::filter(value == TRUE)
    rv
}

#' functions to generate the parts of an upset plot
#'
#' @name generate_X_plot
#' @title generate_X_plot
#' @return a plot of the intersection
#' @param df a long data.frame of data to plot
#' @param id a column in the data for id-ing the individuals
#' @param colour a column to colour data by
#' @param col_vec a named vector of colours for the plot
#' @param intersect_plot plot for the top of the upset
#' @param set_plot plot for the side of the upset
#' @keywords VCF
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate
#' @export
#' @examples
#'
#' generate_intersect_plot(df, id = col)
#' generate_set_plot(df, id = col)
#' generate_group_plot(df, id = col)

generate_intersect_plot <- function(x, ...) {
    UseMethod("generate_intersect_plot")
}

#' @rdname generate_X_plot
#' @export

generate_intersect_plot.data.frame <- function(df,
                                               id,
                                               colour = NULL,
                                               col_vec = NULL,
                                               ...) {
    if ( missing(colour) ) {
        df <- dplyr::group_by(df, {{id}})
    } else {
        df <- dplyr::group_by(df, {{id}}, {{colour}})
    }
    df <- df %>%
        dplyr::summarise(groups = paste0(name, collapse = ", "),
                         .groups = "keep") %>%
        dplyr::group_by(groups) %>%
        dplyr::mutate(n = n())
    if ( missing(colour) ) {
        plot <- ggplot(df, aes(x = reorder(groups, -n)))
    } else {
        plot <- ggplot(df, aes(x = reorder(groups, -n), fill = {{colour}}))
    }
    plot <- plot +
        geom_bar()+
        labs(x = "", y = "intersect size") +
        guides(x = guide_axis(angle = 90)) +
        theme_classic() +
        guides(x = "none", fill = "none") +
        geom_text(aes(y = n + 10,
                      label = n),
                  size = 3)
    if (! is.null(col_vec)) {
        plot <- plot +
            scale_fill_manual(values = col_vec)
    }
    plot
}

#' @rdname generate_X_plot
#' @export

generate_intersect_plot.upsetti_spaghetti <- function(obj, ...) {
    generate_intersect_plot(obj@data, ...)
}

#' @rdname generate_X_plot
#' @export


generate_set_plot <- function(x, ...) {
    UseMethod("generate_set_plot")
}
#' @rdname generate_X_plot
#' @export

generate_set_plot.data.frame <- function(df, colour = NULL, col_vec, ...) {
    df <- df %>%
        dplyr::group_by(name) %>%
        dplyr::mutate(n = n())
    if ( missing(colour) ) {
        plot <- ggplot(df, aes(y = reorder(name, n)))
    } else {
        plot <- ggplot(df, aes(y = reorder(name, n), fill = {{colour}}))
    }
    plot <- plot +
        geom_bar() +
        theme_classic() +
        guides(y = "none", fill = "none") +
        labs(y = "", x = "set size") +
        scale_x_reverse()
    if (! is.null(col_vec)) {
        plot <- plot +
            scale_fill_manual(values = col_vec)
    }
    plot
}

#' @rdname generate_X_plot
#' @export

generate_set_plot.upsetti_spaghetti <- function(obj, ...) {
    generate_set_plot(obj@data, ...)
}

extract_scale <- function(plot, axis, ...) {
    layer_scales(plot)[[axis]]$range$range
}

#' @rdname generate_X_plot
#' @export


generate_group_plot <- function(x, ...) {
    UseMethod("generate_group_plot")
}

#' @rdname generate_X_plot
#' @export

generate_group_plot.data.frame <- function(df,
                                           id,
                                           set_plot,
                                           intersect_plot,
                                           ...) {
    xaxis <- extract_scale(intersect_plot, "x")
    yaxis <- extract_scale(set_plot, "y")
    df %>%
        dplyr::group_by({{id}}) %>%
        dplyr::summarise(groups = paste0(name, collapse = ", ")) %>%
        dplyr::mutate(split = groups) %>%
        tidyr::separate_longer_delim(split, ", ") %>%
        ggplot(aes(
            x = factor(groups, levels = xaxis),
            y = factor(split, levels = yaxis),
            group = groups
        )) +
        geom_point() +
        geom_tile(
            data = expand.grid(x = xaxis, y = yaxis[1:length(yaxis) %% 2 == TRUE]),
            aes(x = x, y = y),
            inherit.aes = FALSE,
            fill = c("grey95")
        ) +
        geom_point(
            data = expand.grid(x = xaxis, y = yaxis),
            aes(x = x, y = y),
            inherit.aes = FALSE,
            colour = "grey80"
        ) +
        geom_point() +
        geom_line() +
        labs(x = "") +
        guides(y = guide_axis(title = ""), x = "none") +
        theme_classic()
}

#' @rdname generate_X_plot
#' @export


generate_group_plot.upsetti_spaghetti <- function(obj, ...) {
    if ("intersect" %in% names(obj@plots) &
        "set" %in% names(obj@plots)) {
        generate_group_plot(obj@data,
                            set_plot = obj@plots$set,
                            intersect_plot = obj@plots$intersect,
                            ...)
    } else {
        generate_group_plot(obj@data, ...)
    }}


#' combine the plot
#'
#' @title process_upsetti_plots
#' @param an upsetti_spaghetti object
#' @return an upsetti_spaghetti object with all 3 plots generated
#' @keywords VCF
#' @import patchwork
#' @export
#' @examples
#'
#' assemble_upset(obj)

process_upsetti_plots <- function(obj, ...) {
    if ( ! "intersect" %in% names(obj@plots) ) {
        obj@plots$intersect <-
            generate_intersect_plot(obj, ...)
    }
    if ( ! "set" %in% names(obj@plots) ) {
        obj@plots$set <-
            generate_set_plot(obj, ...)
    }
    if ( ! "group" %in% names(obj@plots) ) {
        obj@plots$group <-
            generate_group_plot(obj, ...)
    }
    obj
}

#' combine the plot
#'
#' @title assemble_upset
#' @return an upsetti_spaghetti object
#' @keywords VCF
#' @import patchwork
#' @export
#' @examples
#'
#' assemble_upset(obj)

assemble_upset <- function(x, ...) {
    UseMethod("assemble_upset")
}

#' @rdname assemble_upset
#' @export

assemble_upset.upsetti_spaghetti <- function(obj, ...) {
    design = c(
        "
    #AAAAA
    #AAAAA
    BCCCCC
    ")
    obj@plots$final <- obj@plots$intersect +
        obj@plots$set +
        obj@plots$group +
        patchwork::plot_layout(design = design)
    obj
}



#' Take a data.frame and return an upset plot stored in an object of the class upsetti_spaghetti
#'
#' @title plot_upset
#' @return an upsetti_spaghetti object
#' @keywords VCF
#' @export
#' @examples
#'
#' plot_upset(obj)

plot_upset <- function(df, ...) {
    rv <- generate_upsetti(df, ...)
    rv <- process_upsetti_plots(rv, ...)
    rv <- assemble_upset(rv, ...)
    rv
}
