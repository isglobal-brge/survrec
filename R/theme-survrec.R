#' Plot theme used by the survrec graphics
#'
#' A clean `ggplot2` theme shared by all the plotting functions of the
#' package: light horizontal grid, no panel border, muted axis text and
#' room for the legend at the bottom. Exported so that the figures of a
#' report can be restyled consistently. It matches the theme of the
#' companion package `gcmrec`.
#'
#' @param base_size Base font size in points.
#'
#' @return A `ggplot2` theme object, to be added to any plot.
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ group,
#'   data = MMC,
#'   type = "wang-chang"
#' )
#' autoplot(fit) + theme_survrec(base_size = 14)
#' @export
theme_survrec <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.3,
        colour = "grey88"
      ),
      axis.title = ggplot2::element_text(colour = "grey25"),
      axis.text = ggplot2::element_text(colour = "grey35"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.1)
      ),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      plot.caption = ggplot2::element_text(
        colour = "grey45",
        size = ggplot2::rel(0.8)
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(colour = "grey25"),
      strip.text = ggplot2::element_text(face = "bold", colour = "grey25")
    )
}

# Qualitative palette of the package: readable in grayscale and safe for
# the most common types of color blindness (Okabe-Ito). Shared with the
# companion package gcmrec.
survrec_palette <- function(n) {
  pal <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000"
  )
  if (n <= length(pal)) {
    pal[seq_len(n)]
  } else {
    grDevices::colorRampPalette(pal)(n)
  }
}
