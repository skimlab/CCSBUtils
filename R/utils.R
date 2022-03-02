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

#' Combine SummarizedExperiment objects and correct batch effect
#'
#' @param se.list a list of SummarizedExperiment objects
#' @param ID a column name or index for gene ID
#' @param by.assay assay
#' @param ref.batch passed to \code{\link{ComBat}}
#' @param visual.aid if TRUE, add UMAPs before and after batch correction
#' @return combined SummarizedExperiment object
combine_batches <- function(se.list, ID = "ID",
                            by.assay = "data",
                            correct.batch = TRUE,
                            ref.batch = NULL,
                            visual.aid = FALSE,
                            ...) {
  # find common genes first
  common_ids <- rowData(se.list[[1]])[[ID]]
  for (k in 2:length(se.list)) {
    common_ids <- intersect(common_ids, rowData(se.list[[k]])[[ID]])
  }

  # Putting together
  assay_list <- list()

  for (i in seq_along(assays(se.list[[1]]))) {
    an_assay <-
      data.frame(ID = common_ids, row.names = common_ids)

    for (se in se.list) {
      an_assay <-
        left_join(
          an_assay,
          data.frame(ID = rowData(se)[[ID]], assays(se)[[i]]),
          by = "ID"
        )
    }

    rownames(an_assay) <- an_assay[[ID]]
    assay_list[[i]] <- as.matrix(an_assay[-1])
  }
  names(assay_list) <- names(assays(se.list[[1]]))

  clin_data <-
    lapply(se.list,
           FUN = function(se) {
             data.frame(colData(se))
           }
    ) %>%
    bind_rows()

  row_data_df <- data.frame(ID = common_ids)
  rownames(row_data_df) <- row_data_df$ID

  se <- SummarizedExperiment(assays = assay_list)
  rownames(se) <- row_data_df$ID

  rowData(se) <- row_data_df
  colData(se) <- DataFrame(clin_data)

  if (correct.batch) {
    # Quantile normalization
    assay(se, "quantile", withDimnames = FALSE) <-
      normalize.quantiles.robust(assays(se)[[by.assay]])

    # ComBat
    assays(se)[["bc"]] <-
      ComBat(assays(se)[["quantile"]], batch = se$study, ref.batch = ref.batch)
  }

  if (visual.aid) {
    message("visual aid being created....")
    message("\t@metadata$umap_before_batch_correction")
    metadata(se)[["umap_before_batch_correction"]] <-
      compute_umap(se, assay = "quantile") %>% plot_projection(...)
    message("\t@metadata$umap_after_batch_correction")
    metadata(se)[["umap_after_batch_correction"]] <-
      compute_umap(se, assay = "bc") %>% plot_projection(...)
  }

  se
}

#' Combine SummarizedExperiment objects and correct batch effect
#'
#' @param df.list a list of data frames
#' @param ID a column name or index for gene ID
#' @param ref.batch passed to \code{\link{ComBat}}
#' @return combined data frame
combine_batches_dfs <- function(df.list, ID = "ID",
                                ref.batch = NULL,
                                correct.batch = TRUE) {
  # make a named list if no name
  if (is.null(names(df.list))) {
    names(df.list) <- 1:length(df.list)
  }

  # find common genes first
  common_ids <- df.list[[1]][[ID]]
  for (k in 2:length(df.list)) {
    common_ids <- intersect(common_ids, df.list[[k]][[ID]])
  }

  df.combined <-
    lapply(df.list,
           FUN = function(df) {
             df[match(common_ids, df[[ID]]), -match(ID, colnames(df))]
           }) %>%
    bind_cols() %>%
    cbind(ID2 = common_ids, .)

  colnames(df.combined)[1] <- ID

  if (correct.batch) {
    df.combined[-1] <-
      normalize.quantiles.robust(as.matrix(df.combined[-1]))

    batch_names <- names(df.list)
    batch_seq <-
      lapply(
        seq_along(df.list),
        FUN = function(i) {
          rep(batch_names[i], length(df.list[[i]])-1)
        }
      ) %>% unlist()

    # ComBat
    df.combined[-1] <-
      ComBat(df.combined[-1], batch = batch_seq, ref.batch = ref.batch)
  }

  df.combined
}
