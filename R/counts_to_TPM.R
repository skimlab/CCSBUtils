#' A function to extract 'gene_length' from GTF data frame
#'
#' @param gtf_df data.frame formatted GTF file.  Can be obtained by:
#'
#'     gtf <- rtracklayer::import('Homo_sapiens.GRCh38.94.gtf')
#'     gtf_df = as.data.frame(gtf)
#'
#' @return feature_lengths in data frame
gtf_to_gene_length <- function(gtf_df) {
  # internal function
  get_gene_length <- function(df) {
    x <- filter(df, type == "exon")
    min_start <- min(x$start)
    max_end <- max(x$end)
    bps <- rep(0, max_end - min_start+1)
    for (i in 1:nrow(x))
      bps[x$start[i]:x$end[i] - min_start+1] <- 1
    sum(bps)
  }

  gtf_df %>%
    split(.$gene_id) %>%
    map(gtf_to_length) %>%
    stack -> feature_lengths

  feature_lengths %>%
    dplyr::rename(width = values, gene_id = ind)
}

#' Convert counts to transcripts per million (TPM).
#'
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#'
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2
#'
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'
#'    Source for the code: https://gist.github.com/slowkow/c6ab0348747f86e2748b
#'
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength = rep(1, ncol(counts))) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}


#' Convert counts to counts per million (CPM).
#'
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to counts per million.
#'
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @return cpm A numeric matrix normalized by library size
counts_to_cpm <- function(counts) {
  cpm = counts/colSums(counts)

  # Copy the row and column names from the original matrix.
  colnames(cpm) <- colnames(counts)
  rownames(cpm) <- rownames(counts)
  return(cpm)
}
