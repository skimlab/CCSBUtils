
#' A function to perform *UMAP* projection of data
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param assay A name of the slot data is stored
#' @return A \code{\link{SummarizedExperiment}}  object with *UMAP* projection stored in \code{metadata(se[[umap]]} slot
#' @export
compute_umap <- function(se, assay = "data") {
  S4Vectors::metadata(se)[["umap"]] <- data.frame(uwot::umap(t(SummarizedExperiment::assays(se)[[assay]])),
                                                  row.names = colnames(se))
  colnames(S4Vectors::metadata(se)[["umap"]]) <- c("umap_1", "umap_2")
  se
}

#' A function to perform *t-SNE* projection of data
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param assay A name of the slot data is stored
#' @return A \code{\link{SummarizedExperiment}} object with *t-SNE* projection stored in \code{metadata(se[[tsne]]} slot
#' @export
compute_tsne <- function(se, assay = "data", num_threads = 1) {
  S4Vectors::metadata(se)[["tsne"]] <- data.frame(Rtsne::Rtsne(t(SummarizedExperiment::assays(se)[[assay]]), num_threads = num_threads)$Y,
                                                  row.names = colnames(se))
  colnames(S4Vectors::metadata(se)[["tsne"]]) <- c("tsne_1", "tsne_2")
  se
}

#' A function to perform *PCA* projection of data
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param assay A name of the slot data is stored
#' @return A \code{\link{SummarizedExperiment}} object with *PCA* projection stored in \code{metadata(se[[pca]]} slot
#' @export
compute_pca <- function(se, assay = "data") {
  S4Vectors::metadata(se)[["pca"]] <- data.frame(stats::prcomp(t(SummarizedExperiment::assays(se)[[assay]]))$x,
                                                 row.names = colnames(se))
  colnames(S4Vectors::metadata(se)[["pca"]]) <- c("pc_1", "pc_2")
  se
}

#' A function to visualize/plot a projection
#'
#' @param se A \code{\link{SummarizedExperiment}} object
#' @param reduction A name of the slot data where projection is is stored: 'umap', 'tsne' or 'pca'
#' @param feature_color If specified, the projection will be color-coded using the specified feature from \code{colData(se)}
#' @param feature_shape If specified, the projection will be shape-coded using the specified feature from \code{colData(se)}
#' @return A \code{\link{ggplot}} object
#' @export
plot_projection <- function(se, reduction = "umap", feature_color = NA, feature_shape = NA) {
  reduction <- S4Vectors::metadata(se)[[reduction]]
  projection_df <- data.frame(name.xy = rownames(SummarizedExperiment::colData(se)),
                              x = reduction[[1]],
                              y = reduction[[2]])
  # to keep all other column data and make them available
  projection_df <- cbind(
    projection_df,
    data.frame(colData(se))
  )
  if (!is.na(feature_color)) {
    if (!is.na(feature_shape)) {
      projection_df <- cbind(projection_df, feature_color = se[[feature_color]], feature_shape = se[[feature_shape]])
      gp <-
        ggplot2::ggplot(projection_df, ggplot2::aes(x = x, y = y,
                                                    color = feature_color, shape = feature_shape, label = name.xy))
    } else {
      projection_df <- cbind(projection_df, feature_color = se[[feature_color]])
      gp <-
        ggplot2::ggplot(projection_df, ggplot2::aes(x = x, y = y,
                                                    color = feature_color, label = name.xy))
    }
  } else {
    gp <-
      ggplot2::ggplot(projection_df, ggplot2::aes(x = x, y = y, label = name.xy))
  }

  gp +
    ggplot2::geom_point() +
    ggplot2::xlab(colnames(reduction)[1]) +
    ggplot2::ylab(colnames(reduction)[2]) +
    ggplot2::labs(color = feature_color, shape = feature_shape)
}
