#' predict outcomes of observations in data using *selected features*
#'
#' @param model a CV trained model, an output of cv_loop_train_iter
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an observation vector.
#'             Since this is an input to \code{\link{glmnet}}, it should be the format that can be used
#'             with \code{\link{glmnet}}
#' @param type 'raw' (default) or 'prob'
#' @return Either outputs class labels ('raw) or
#'         probability tables, a column for each class label ('prob')
predict.type <- function(model, data, type = "raw") {
  predict(model$fit, select(data, model$selected_features$feature), type = type)
}


#' predict outcomes of observations in data using *selected features*
#'
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an observation vector.
#'             Since this is an input to \code{\link{glmnet}}, it should be the format that can be used
#'             with \code{\link{glmnet}}
#' @param cls class labels
#' @param model a CV trained model, an output of cv_loop_train_iter
#' @return A list of
#'     class (class labels),
#'     trained (names of trained samples),
#'     row.names (selected features),
#'     pred (predicted outcomes), and
#'     prob (probabilities of predicted outcomes)
prediction_tables <- function(data, cls, model) {
  ret <- list()
  ret[["class"]] <- cls
  ret[["trained"]] <- rownames(data) %in% models$train_samples
  ret[["row.names"]] <- rownames(data)
  ret[["pred"]] <-
    data.frame(pred = predict.type(model, data),
               prob = predict.type(model, data, type = "prob"))
  ret
}

#' Generate a list of prediction summary
#'
#' @param pred_tbls An output from \code{\link{prediction_tables}}
#' @return A list of prediction summary
prediction_summary <- function(pred_tbls) {

  ret <- list(class = pred_tbls$class, trained = pred_tbls$trained, row.names = pred_tbls$row.names,
              summary = list())

  ret$summary$p.table <- pred_tbls$pred
  ret$summary$confusion <- confusionMatrix(pred_tbls$pred$pred, pred_tbls$class)
  ret$summary$roc <- pROC::roc(pred_tbls$class, pred_tbls$pred[[2]])  # prob.[type]

  ret$summary$accuracy <- ret$summary$confusion$overall[["Accuracy"]]
  ret$summary$kappa <- ret$summary$confusion$overall[["Kappa"]]
  ret$summary$auc <- auc(ret$summary$roc)
  ret$summary$F1 <- ret$summary$confusion$byClass[["F1"]]
  ret$summary$NPV <- ret$summary$confusion$byClass[["Neg Pred Value"]]
  ret$summary$PPV <- ret$summary$confusion$byClass[["Pos Pred Value"]]
  ret$summary$sensitivity <- ret$summary$confusion$byClass[["Sensitivity"]]
  ret$summary$specificity <- ret$summary$confusion$byClass[["Specificity"]]

  ret
}

#' Generate a summary of performance of CV trained models
#'
#' @param cv_trained An output of \code{\link{cv_loop_train}} or \code{\link{cv_loop_train_parallel}}
models_performance_summary <-
  function(cv_trained) {

    data <- cv_trained$data
    cls <- cv_trained$class
    models <- cv_trained$models

    K <- length(models)

    train_ <- list()
    test_ <- list()
    for (k in 1:K) {
      data_train <- data[models[[k]]$train_index, ]
      data_test <- data[-(models[[k]]$train_index), ]

      class_train <- cls[models[[k]]$train_index]
      class_test <- cls[-(models[[k]]$train_index)]

      prediction_tables(data = data_train, cls = class_train, model = models[[k]]) %>%
        prediction_summary() -> train_[[k]]

      prediction_tables(data = data_test, cls = class_test, model = models[[k]]) %>%
        prediction_summary() -> test_[[k]]
    }

    list(train = train_, test = test_)
  }

#' Extract performance metrics from \code{performance_summary} slot
#'
#' @param p.summary Either "train" or "test" slot of an output of \code{\link{models_performance_summary}}
#' @param slot An index
extract_metrics <- function(p.summary) {
  data.frame(
    idx = 1:length(p.summary),
    accuracy = sapply(p.summary, function(x) x$summary[["accuracy"]]),
    kappa = sapply(p.summary, function(x) x$summary[["kappa"]]),
    roc = sapply(p.summary, function(x) x$summary[["roc"]]),
    auc = sapply(p.summary, function(x) x$summary[["auc"]]),
    F1 = sapply(p.summary, function(x) x$summary[["F1"]]),
    NPV = sapply(p.summary, function(x) x$summary[["NPV"]]),
    PPV = sapply(p.summary, function(x) x$summary[["PPV"]]),
    sensitivity = sapply(p.summary, function(x) x$summary[["sensitivity"]]),
    specificity = sapply(p.summary, function(x) x$summary[["specificity"]])
  )
}

#' Extract all performance metrics
#'
#' @param p.summary Either "train" or "test" slot of an output of \code{\link{models_performance_summary}}
extract_metrics_all <- function(p.summary) {
  # p.summary is a array (K trials) of performance summary.
  # nn gets if both "uniform" and "weighted" are tried.
  nn <- names(p.summary[[1]]$summary)
  r <- list()
  for (n in nn) {
    r[[n]] <- extract_metrics(p.summary, slot = n)
  }

  bind_rows(r, .id = "feature_weight")
}


#' Workflow to get all classification performance
#'
#' @param cv_trained An output of \code{\link{cv_loop_train}} or \code{\link{cv_loop_train_parallel}}
#' @return A list of classification performance summary
classification_summary_workflow <- function(cv_trained, show_feature_map = FALSE) {
  cv_trained_summary <- list(cv_trained = cv_trained)

  cv_trained_summary[["performance_summary"]] <- models_performance_summary(cv_trained)
  cv_trained_summary[["accuracy"]] <-
    bind_rows(
      cbind(feature_weight = cv_trained$models[[1]]$feature_weights,
            type = "train",
            extract_metrics(cv_trained_summary[["performance_summary"]][["train"]])),
      cbind(feature_weight = cv_trained$models[[1]]$feature_weights,
            type = "test",
            extract_metrics(cv_trained_summary[["performance_summary"]][["test"]]))
    )

  # print(cv_trained_summary[["accuracy"]])

  cv_trained_summary[["accuracy"]] %>%
    pivot_longer(-c("type", "feature_weight", "idx")) %>%
    ggplot(aes(name, value)) +
    geom_jitter(aes(color = name)) +
    xlab("performance metrics") ->  accuracy_gp


  wrap_up_accuracy_gp <- function(gp) {
    gp  +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      theme(legend.position = "none") +
      facet_wrap(c("feature_weight", "type"), nrow = 1)
  }

  # violin plots
  wrap_up_accuracy_gp (
    accuracy_gp +
      geom_violin(aes(color = name))
  ) -> cv_trained_summary[["accuracy_violin"]]

  # box plots
  wrap_up_accuracy_gp (
    accuracy_gp +
      geom_boxplot(aes(color = name), notch = TRUE, width = 0.8)
  ) -> cv_trained_summary[["accuracy_boxplot"]]

  # print(cv_trained_summary[["accuracy_boxplot"]])

  # summary features
  cv_trained_summary[["feature_maps"]] <- construct_feature_maps(cv_trained$models)

  cv_trained_summary <- construct_feature_heatmaps(cv_trained_summary)

  # temp variable to ease the coding below
  feature_maps <- cv_trained_summary[["feature_maps"]]

  # count the number of features/predictors
  n_features <- function(x_map) {
    colSums(select(x_map, starts_with("X")) > 0)
  }

  n_features_list <- list()
  if (!is.null(feature_maps[["uniform"]])) {
    n_features_list[["uniform.n_features"]] <- n_features(feature_maps[["uniform"]]$feature_map)
    n_features_list[["uniform.n_predictors"]] <- n_features(feature_maps[["uniform"]]$predictor_map)

    if (show_feature_map)
      print(feature_maps[["uniform"]]$coef_heatmap)

  }

  if (!is.null(feature_maps[["weighted"]])) {
    n_features_list[["weighted.n_features"]] <- n_features(feature_maps[["weighted"]]$feature_map)
    n_features_list[["weighted.n_predictors"]] <- n_features(feature_maps[["weighted"]]$predictor_map)

    if (show_feature_map)
      print(feature_maps[["weighted"]]$coef_heatmap)
  }

  cv_trained_summary[["n_features"]] <- n_features_list

  cv_trained_summary
}


#' Generate a consensus summary of outcomes
#'
#' @param cv_trained An output of \code{\link{cv_loop_train}} or \code{\link{cv_loop_train_parallel}}
#' @return A table of consensus summary of outcomes
prediction_consensus_summary <- function(cv_trained) {
  lapply(cv_trained$models,
         function(m) {
           predict.type(m, cv_trained$data)
         }) -> x
  names(x) <- sprintf("X%03d", 1:length(cv_trained$models))

  lapply(cv_trained$models,
         function(m) {
           bind_cols(sample_id = rownames(cv_trained$data),
                     predict.type(m, cv_trained$data, type = "prob"))
         }) -> x.prob
  names(x.prob) <- sprintf("X%03d", 1:length(cv_trained$models))

  x2 <- bind_cols(ID = rownames(cv_trained$data),
                  cls = cv_trained$class,
                  x)

  cls_labels = levels(cv_trained$class)
  for (cls in cls_labels) {
    x2[[cls]] <- rowSums(x2[, -c(1:2)] == cls)
  }

  x2[["majority"]] <- cls_labels[apply(x2[cls_labels], MARGIN = 1, which.max)]
  x2[["purity"]] <- apply(x2[cls_labels], MARGIN = 1, max)/length(cv_trained$models)

  x2 <- mutate(x2, mistake = cls != majority)

  x2 %>%
    #    select(c(1:2), cls_labels, majority, purity, mistake, starts_with("X")) %>%
    select(-starts_with("X"), starts_with("X")) %>%
    arrange(-mistake, purity)
}
