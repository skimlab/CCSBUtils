
#' A function to perform *UMAP* projection of data
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param assay A name of the slot data is stored
#' @return A \code{\link{SummarizedExperiment}}  object with *UMAP* projection stored in \code{metadata(se[[umap]]} slot
#' @export
compute_umap <- function(se, assay = "data") {
  S4Vectors::metadata(se)[["umap"]] <- data.frame(uwot::umap(t(SummarizedExperiment::assays(se)[[assay]])))
  se
}

#' A function to perform *t-SNE* projection of data
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param assay A name of the slot data is stored
#' @return A \code{\link{SummarizedExperiment}} object with *t-SNE* projection stored in \code{metadata(se[[tsne]]} slot
#' @export
compute_tsne <- function(se, assay = "data") {
  S4Vectors::metadata(se)[["tsne"]] <- data.frame(Rtsne::Rtsne(t(SummarizedExperiment::assays(se)[[assay]]))$Y)
  se
}

#' A function to perform *PCA* projection of data
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param assay A name of the slot data is stored
#' @return A \code{\link{SummarizedExperiment}} object with *PCA* projection stored in \code{metadata(se[[pca]]} slot
#' @export
compute_pca <- function(se, assay = "data") {
  S4Vectors::metadata(se)[["pca"]] <- data.frame(stats::prcomp(t(SummarizedExperiment::assays(se)[[assay]]))$x)
  se
}

#' A function to visualize/plot a projection
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param reduction A name of the slot data where projection is is stored: 'umap', 'tsne' or 'pca'
#' @param feature If specified, the projection will be color-coded using the specified feature from \code{colData(se)}
#' @return A \code{\link{ggplot}} object
#' @export
plot_projection <- function(se, reduction = "umap", feature = NA) {
  projection_df <- data.frame(name = rownames(SummarizedExperiment::colData(se)),
                              x = S4Vectors::metadata(se)[[reduction]][[1]],
                              y = S4Vectors::metadata(se)[[reduction]][[2]])
  # to keep all other column data and make them available
  projection_df <- cbind(
    projection_df,
    colData(se)
  )
  if (!is.na(feature)) {
    projection_df <- cbind(projection_df, feature = se[[feature]])
    gp <-
      ggplot2::ggplot(projection_df, ggplot2::aes(x = x, y = y,
                                                  color = feature, label = name))
  } else {
    gp <-
      ggplot2::ggplot(projection_df, ggplot2::aes(x = x, y = y, label = name))
  }

  gp +
    ggplot2::geom_point() +
    ggplot2::xlab(paste(reduction, "1", sep = "_")) +
    ggplot2::ylab(paste(reduction, "2", sep = "_")) +
    ggplot2::guides(color=ggplot2::guide_legend(title=feature))
}
