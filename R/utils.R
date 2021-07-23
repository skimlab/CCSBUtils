#' Split SummarizedExperiment object
#'
#' @param se SummarizedExperiment
#' @param col integer (index) or string (column heading) of colData(se)
#' @return a list of split objects
split_se <- function(se, col) {
  cx <- unique(se[[col]])

  res <- lapply(cx, function(cc) {
    idx <- which(se[[col]] == cc)
    si <- se[, idx]
    metadata(si) <- subset_metadata(se, idx)
    si
  })

  names(res) <- cx
  res
}

#' Subset metadata
#' @noRd
subset_metadata <- function(se, idx) {
  mse <- metadata(se)
  res <- lapply(mse,
                function(df) {
                  if ("data.frame" %in% class(df))
                    df[idx,]
                  else
                    df
                })
  names(res) <- names(mse)
  res
}
