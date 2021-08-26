#' An internal function to construct heatmaps from feature maps
#'
#' @param cv_trained_summary An intermediate struct that contains
#'   \code{feature_maps} and \code{accuracy} slots
construct_feature_heatmaps <- function(cv_trained_summary) {
  coef_pheatmap <- function(m, ...) {
    m <- as.matrix(m)
    ComplexHeatmap::pheatmap(
      m,
      name = "score",
      show_colnames = F,
      show_rownames = F,
      ...
    )
  }

  feature_pheatmap <- function(m, ...) {
    m <- as.matrix(m)
    ComplexHeatmap::pheatmap(
      m,
      name = "score",
      show_colnames = F,
      show_rownames = F,
      ...
    )
  }

  coef_heatmap <- function(m, ...) {
    m <- as.matrix(m)
    mx <- max(abs(m)) * .75
    Heatmap(
      m,
      name = "coef",
      show_column_names = F,
      show_row_names = F,
      col = circlize::colorRamp2(
        breaks = seq(-mx, mx, length.out = 7),
        colors = brewer.pal(n = 7, name = "RdYlBu")
      ),
      # col = rev(brewer.pal(n = 7, name = "RdYlBu")),
      ...
    )
    # ComplexHeatmap::pheatmap(m, name = "score", show_colnames = F, show_rownames = F, ...)
  }

  feature_heatmap <- function(m, ...) {
    m <- as.matrix(m)
    Heatmap(
      m,
      name = "score",
      show_column_names = F,
      show_row_names = F,
      # col = circlize::colorRamp2(breaks = c(0,  1), colors = c("white", "navy")),
      col = rev(brewer.pal(n = 7, name = "RdYlBu")),
      ...
    )
    # ComplexHeatmap::pheatmap(m, name = "score", show_colnames = F, show_rownames = F, ...)
  }

  feature_maps <- cv_trained_summary[["feature_maps"]]

  feature_rowAnnot <-
    rowAnnotation(
      feature = feature_maps$feature_map$overall_score,
      col = list(feature = circlize::colorRamp2(
        breaks = c(0, 1),
        colors = c("white", "darkred")
      ))
    )
  predictor_rowAnnot <-
    rowAnnotation(
      predictor = anno_barplot(feature_maps$predictor_map$overall_score),
      feature = feature_maps$predictor_map$overall_feature_score,
      col = list(
        predictor = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "darkgreen")
        ),
        feature = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "darkred")
        )
      )
    )
  coef_rowAnnot <-
    rowAnnotation(
      coef = anno_barplot(feature_maps$coef_map[["overall_score"]],
                          axis_param = list(side = "top")),
      coef_mad = anno_barplot(feature_maps$coef_map[["overall_score_mad"]],
                              axis_param = list(side = "top")),
      feature = feature_maps$coef_map[["overall_feature_score"]],
      miRNA = anno_text(feature_maps$coef_map$feature),
      col = list(
        coef = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "darkgreen")
        ),
        coef2 = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "darkgreen")
        ),
        feature = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "darkred")
        )
      )
    )

  cv_trained_summary[["accuracy"]] %>%
    filter(feature_weight == "uniform", type == "test") %>%
    dplyr::select(accuracy) -> acc

  feature_colAnnot <-
    HeatmapAnnotation(n_features = anno_barplot(colSums(
      dplyr::select(feature_maps$feature_map, starts_with("X")) > 0
    ),
    axis_param = list(side = "left")))
  predictor_colAnnot <-
    HeatmapAnnotation(
      accuracy = anno_barplot(acc, col = "blue"),
      n_predictors = anno_barplot(colSums(
        dplyr::select(feature_maps$predictor_map, starts_with("X")) > 0
      ),
      axis_param = list(side = "left"))
    )

  feature_maps[["feature_heatmap"]] <-
    feature_heatmap(
      dplyr::select(feature_maps$feature_map, starts_with("X")),
      right_annotation = feature_rowAnnot,
      top_annotation = feature_colAnnot,
      row_order = order(-feature_maps$feature_map$overall_score)
    )

  feature_maps[["predictor_heatmap"]] <-
    feature_heatmap(
      dplyr::select(feature_maps$predictor_map, starts_with("X")),
      right_annotation = predictor_rowAnnot,
      top_annotation = predictor_colAnnot,
      row_order = order(-feature_maps$predictor_map$overall_score)
    )

  feature_maps[["coef_heatmap"]] <-
    coef_heatmap(
      sign(dplyr::select(
        feature_maps$coef_map, starts_with("X")
      )),
      right_annotation = coef_rowAnnot,
      top_annotation = predictor_colAnnot,
      row_order = order(-feature_maps$coef_map$overall_score)
    )

  # feature_maps[["predictor_heatmap2"]] <-
  #   feature_pheatmap(
  #     select(feature_maps$coef_map, starts_with("X")),
  #     right_annotation = predictor_rowAnnot,
  #     top_annotation = predictor_colAnnot,
  #     row_order = order(-feature_maps$predictor_map$overall_score)
  #   )

  # this heatmap requires a bit more work
  x <- dplyr::select(feature_maps$coef_map, starts_with("X"))
  x <-
    t(scale(
      t(x),
      center = F,
      scale = apply(abs(x), 1, function(ax)
        max(ax) * 0.75)
    ))
  feature_maps[["coef_heatmap2"]] <-
    coef_heatmap(
      x,
      right_annotation = coef_rowAnnot,
      top_annotation = predictor_colAnnot,
      row_order = order(-feature_maps$coef_map$overall_score)
    )

  cv_trained_summary[["feature_maps"]] <- feature_maps
  cv_trained_summary
}


#' A function to construct boxplots for selected features
#'
#' @param dd A \code{matrix} or \code{data.frame}
#'           where each row is an observation
#' @param cls A set of classes, either in \code{factor}
#' @param df fold difference, computed by \code{compute_difference}
#' @param order.by which feature to order boxplots, typically `mean` or `diff`
feature_boxplots <-
  function(dd,
           cls,
           df = data.frame(),
           order.by = "") {
    dplyr::bind_cols(cls = cls, dd = data.frame(dd, check.names = F)) %>%
      tidyr::pivot_longer(-cls) ->
      data_cls_long

    data_cls_long %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(median = median(value), mean = mean(value)) %>%
      dplyr::right_join(data_cls_long, by = c("name")) -> data_cls_long

    data_cls_long %>%
      dplyr::left_join(df, by = c('name' = colnames(df)[1])) ->
      data_cls_long

    purrr::map(colnames(df[-1]), function(feature) {
      data_cls_long %>%
        dplyr::mutate(name = fct_reorder(name, eval(parse(text = order.by)))) %>%
        dplyr::select(name, feature) %>%
        dplyr::distinct() %>%
        ggplot(aes(x = name, y = eval(parse(text = feature)))) +
        geom_bar(stat = "identity") +
        ylab(feature) +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    }) -> gp_list
    names(gp_list) <- colnames(df[-1])

    data_cls_long %>%
      dplyr::mutate(name = fct_reorder(name, eval(parse(text = order.by)))) %>%
      ggplot(aes(x = name, y = value, fill = cls)) +
      geom_boxplot(notch = T) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) -> gp_boxplots
    gp_list[["boxplots"]] <- gp_boxplots

    ggarrange(
      plotlist = gp_list,
      ncol = 1,
      nrow = length(gp_list),
      align = "v",
      heights = c(rep(1, length(gp_list) - 1), 8),
      legend = "bottom"
    )
  }


#' Create a heatmap of data with predicted outcomes at the top and consensus
#'
#' @param se \code{\link{SummarizedExperiment}} container
#' @param assay a name of assay slot in \code{se}
#' @param cls_results An output of \code{\link{classification_summary_workflow}}
#' @param scale If TRUE, the data (=assays(se)[[assay]]) will be scaled
#' @param ...  Additional parameters to be passed to \code{\link{Heatmap}}
#' @return Heatmap object
heatmap_prediction <-
  function(se,
           assay = "data",
           cls_results,
           scale = FALSE,
           ...) {
    cls_results[["consensus"]] %>%
      dplyr::select(ID, cls, majority, purity) %>%
      mutate(ord =  purity * (2 * (as.integer(cls) - 1) - 1)) %>%
      as.data.frame ->
      annot_col_df
    rownames(annot_col_df) <- annot_col_df$ID

    n_cls <- length(levels(annot_col_df$cls))
    cls_color <-
      brewer.pal(ifelse(n_cls < 3, 3, n_cls), "Set1")[1:n_cls]
    names(cls_color) <- levels(annot_col_df$cls)

    top_annot_color <-
      list(cls = cls_color,
           majority = cls_color,
           purity = colorRamp2(c(0.5, 1.0), c("white", "black")))

    dd <- assays(se)[[assay]]
    dd <- dd[, annot_col_df$ID]

    annot_col <-
      HeatmapAnnotation(df = annot_col_df[c("cls", "majority", "purity")],
                        col = top_annot_color[c("cls", "majority", "purity")],
                        which = "column")

    if (scale) {
      dd <- t(scale(t(dd)))
    }

    # making sure column names retained
    colnames(dd) <- annot_col_df$ID

    annot_row <-
      HeatmapAnnotation(
        diff = compute_difference(t(dd), annot_col_df$cls),
        col = list(diff = colorRamp2(c(-2, 0, 2), c(
          "blue", "white", "red"
        ))),
        which = "row"
      )

    ComplexHeatmap::Heatmap(
      as.matrix(dd),
      name = assay,
      top_annotation = annot_col,
      right_annotation = annot_row,
      column_split = annot_col_df$cls,
      column_order = order(annot_col_df$ord),
      ...
    )
  }


#' Create a heatmap of projected data (coef x data) with predicted outcomes at
#' the top and consensus
#'
#' @param se \code{\link{SummarizedExperiment}} container
#' @param assay a name of assay slot in \code{se}
#' @param cls_results An output of \code{\link{classification_summary_workflow}}
#' @param cv_model output of caret::train or glmnet::cv_glmnet
#' @param scale If TRUE, the data (=assays(se)[[assay]]) will be scaled
#' @param ...  Additional parameters to be passed to \code{\link{Heatmap}}
#' @return Heatmap object
heatmap_glmnet <-
  function(se,
           assay = "data",
           cls_results,
           cv_model,
           scale = FALSE,
           ...) {
    cls_color <- c("lightgray", "black")
    names(cls_color) <- unique(cls_results$cv_trained$class)

    top_annot_color <-
      list(majority = cls_color,
           cls = cls_color)

    purity_colorRamp <- colorRamp2(c(0.5, 1.0), c("white", "black"))
    top_annot_color[["purity"]] <- purity_colorRamp


    cls_results[["consensus"]] %>% dplyr::select(ID, cls, majority, purity) %>%
      mutate(ord =  purity * (2 * (as.integer(cls) - 1) - 1)) %>%
      as.data.frame ->
      annot_top_df
    # rownames(annot_top_df) <- annot_top_df$sample_uid
    rownames(annot_top_df) <- annot_top_df$ID


    annot_top <-
      HeatmapAnnotation(df = annot_top_df[c("cls", "majority", "purity")],
                        col = top_annot_color[c("cls", "majority", "purity")],
                        which = "column")

    which.beta <-
      which(cv_model$finalModel$lambda == cv_model$finalModel$lambdaOpt)

    bb <- cv_model$finalModel$beta[, which.beta]
    a0 <- cv_model$finalModel$a0[which.beta]

    dd <- assays(se)[[assay]]
    dd <- dd[names(bb), annot_top_df$ID]

    prob <- predict(cv_model, t(dd), type = "prob")
    prob <- prob[[2]]
    # prob <- plogis(colSums(dd*bb)+a0)


    annot_bottom <-
      columnAnnotation(prob = anno_points(prob, gp = gpar(col = ifelse(
        prob > 0.5, "red", "blue"
      ))))

    annot_right <-
      rowAnnotation(coef = bb, col = list(coef = colorRamp2(
        c(-max(abs(bb)), 0, max(abs(bb))), c("blue", "white", "red")
      )))

    dx <- dd * bb + a0
    if (scale) {
      dx <- t(scale(t(dx)))
    }

    ComplexHeatmap::Heatmap(
      dx,
      name = "projected",
      top_annotation = annot_top,
      bottom_annotation = annot_bottom,
      right_annotation = annot_right,
      column_split = ifelse(colSums(dd * bb) + a0 > 0, 'pos', 'neg'),
      column_order = order(prob),
      row_order = order(bb),
      ...
    )
  }


#' Create line plots to show train trajectories
#'
#' @param cv_model an output of caret::train or glmnet::cv_glmnet
#' @return \code{\link{ggplot2}} object
ggplot_train_trajectory <- function(cv_model) {
  dplyr::left_join(cv_model$results,
                   data.frame(cv_model$finalModel[c("lambda", "df")]),
                   by = "lambda") %>%
    tidyr::pivot_longer(c("Accuracy", "df")) %>%
    ggplot2::ggplot(aes(x = lambda, y = value)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::xlab("Regularization parameter") +
    ggplot2::ylab("") +
    ggplot2::facet_grid(name ~ ., scales = "free")
}

#' Create a heatmap of correlation between features (rows)
#'
#' @param dd data (row: feature, column: observation)
#' @param importance a vector of `importance` values of features
#' @return \code{\link{ComplexHeatmap}} object
heatmap_features_cor <- function(dd, importance, ...) {
  dd %>% t %>% cor %>%
    Heatmap(
      right_annotation = rowAnnotation(importance = anno_barplot(importance,
                                                                 axis_param = list(side = "top"))),
      name = "correlation",
      col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "yellow", "red")),
      ...
    )
}
