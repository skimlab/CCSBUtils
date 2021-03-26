width.Heatmap <- function(ht, unitTo = "inch") {
  grid::convertWidth(
    sum(component_width(ht)) + unit(4, units = "mm"),
    unitTo = unitTo,
    valueOnly = TRUE
  )
}

height.Heatmap <- function(ht, unitTo = "inch") {
  grid::convertHeight(
    sum(component_height(ht)) + unit(4, units = "mm"),
    unitTo = unitTo,
    valueOnly = TRUE
  )
}

#' A function to estimate the drawing size of Heatmap object
#'
#' @param ht Heatmap object, an output of draw(a_Heatmap_object)
#'        and \code{accuracy} slots
#' @param unitTo The coordinate system to convert the unit to. See the unit function for valid coordinate systems,
#'               to be passed to grid::convertHeight() and grid::convertWidth() functions.
#' @return a list of width and height in "unitTo"
size.Heatmap <- function(ht, unitTo = "inch") {
  list(width = width.Heatmap(ht, unitTo),
       height = height.Heatmap(ht, unitTo))
}

#' A function to save a Heatmap object to a PDF file
#'
save.PDF.Heatmap <- function(ht,
                             filename,
                             width = width.Heatmap(ht, "inch"),
                             height = height.Heatmap(ht, "inch")) {
  pdf(file = filename, width = width, height = height)
  draw(ht)
  dev.off()
}

save.PNG.Heatmap <- function(ht,
                             filename,
                             width = width.Heatmap(ht, "inch"),
                             height = height.Heatmap(ht, "inch"),
                             units = "in",
                             res = 100) {
  png(file = filename, width = width, height = height,
      units = units, res = res)
  draw(ht)
  dev.off()
}
