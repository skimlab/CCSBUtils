#' Compute Robust Z score
#'
#' @param dd a matrix
#' @return robust z score in matrix
robust_zscore <- function(dd) {
  require(dplyr)

  dd_median <- apply(dd, MARGIN=1, FUN=median)
  dd_mad <- apply(dd, MARGIN=1, FUN=mad)

  rz <- sweep(dd, 1, dd_median, "-")
  rz <- sweep(rz, 1, dd_mad, "/")

  rz
}

#' Quantize values in a matrix
#'
#' @param rz robust z-score matrix
#' @return quantized robust z-score matrix
quantize_robust_zscore <- function(rz) {
  # quantization,
  #    1 if larger than median + 1 MAD,
  #   -1 if smaller than median - 1 MAD
  #    0 otherwise
  #
  rz <- as.matrix(rz)

  rz[is.nan(rz)] <- 0
  rz[is.infinite(rz)] <- sign(rz[is.infinite(rz)])
  rz[rz > 1] <-  1
  rz[rz < -1] <- -1
  rz[rz > -1 & rz < 1] <- 0

  rz
}
