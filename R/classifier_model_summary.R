#' predict outcomes of observations in data using *selected features*
#'
#' @param model a CV trained model, an output of cv_loop_train_iter
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param type 'raw' (default) or 'prob'
#' @return a named vector, either outputs class labels ('raw) or probability
#'   tables, a column for each class label ('prob')
predict.type <- function(model, data, type = "raw") {
  features <- colnames(model$fit$trainingData)
  features <- features[-length(features)]

  if (sum(is.na(match(features, colnames(data)))) > 0)
    stop("some features are missing in the data...")

  pred <- predict(model$fit, data[, features], type = type)
  if (type == "raw")
    names(pred) <- rownames(data)
  else {
    pred <- data.frame(Sample_ID = rownames(data), pred)
    pred[["prob"]] <- pred[[model$fit$levels[1]]]
    pred[["predicted"]] <-
      predict(model$fit, data[, features], type = "raw")
    rownames(pred) <- rownames(data)
  }

  pred
}


#' apply a list of models to a set of samples
#'
#' returns a table where each rows is predicted outcomes of a sample, by
#' applying predictive models
#'
#' @param models a list of models which is an outcome of
#'   \code{\link{cv_loop_train_iter}} or \code{\link{cv_train_final}}, a model
#'   trained for a class
#' @param dd   input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @returns a table where each rows is predicted outcomes of applying predictive
#'   models to a sample
predict.consensus <- function(models, data) {
  lapply(1:length(models),
         function(k)
           predict.type(models[[k]], data = data, type = "raw")) %>%
    bind_rows %>%
    lapply(table) %>%
    bind_rows(.id = "Sample_ID") %>%
    data.frame -> pred

  cpred <- colnames(pred)[-1]
  pred[["prob"]] <- pred[[cpred[1]]] / length(models)
  pred[["predicted"]] <-
    cpred[apply(pred[, cpred], MARGIN = 1, which.max)]

  rownames(pred) = pred$Sample_ID
  pred
}


#
# OBSOLETE, no use
#
#' convert the outcomes of \code{\link{predict.consensus}}
#'   to the one equivalent of \code{\link{predict.type}}
#'
#' @param cpred an outcome of \code{\link{predict.consensus}}
#' @return the last column is added with the multi-class labels (flatted)
predict.consensus.type <- function(cpred) {
  res <- cpred$predicted
  names(res) <- cpred$Sample_ID
  res
}


#' get subset of samples with 'predicted' = 'cls'
#'
#' @param pred a table with 'Sample_ID' and 'predicted' columns
#' @param cls target class label
#' @return a list of 'Sample_ID's
get.subsamples <- function(pred, cls) {
  # pred %>%
  #   filter(predicted == cls) -> x
  # x$Sample_ID
  filter(pred, predicted == cls)[["Sample_ID"]]
}



#
# OBSOLETE, no use
#
#' apply a list of models to a set of samples
#'
#' returns a table where each rows is predicted outcomes of a sample, by
#' applying predictive models
#'
#' @param models a list of models which is an outcome of
#'   \code{\link{cv_loop_train_iter}} or \code{\link{cv_train_final}}, a model
#'   trained for a class
#' @param dd   input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param type either "raw" or "prob", for the number/class predictions or class
#'   probabilities, respectively. Class probabilities are not available for all
#'   classification models
#' @returns a table where each rows is predicted outcomes of applying predictive
#'   models to a sample
predict.multiClass <- function(models, dd, type = "raw") {
  # models is a list of model: fit, selected_features
  if (is.null(names(models)))
    names(models) <- 1:length(models)

  predicted <-
    lapply(models,
           function(m, d) {
             features <- colnames(m$fit$trainingData)
             features <- features[-length(features)]

             predict(m$fit, dplyr::select(d, feature), type = type)
           },
           dd)
  tibble(ID = rownames(dd), do.call(bind_cols, predicted)[names(models)])
}


#' apply a list of models to a set of samples
#'
#' returns a table where each rows is predicted outcomes of a sample, by
#' applying predictive models
#'
#' @param models a list of models each of which is an outcome of
#'   \code{\link{cv_loop_train}}
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param final.class  the name (column heading) of final class outcomes,
#'   majority among multiple classes
#' @param min.prob minimum probability for a sample to be classified to a class
#'   (default: 0.25)
#' @returns a table where each rows is predicted outcomes of applying predictive
#'   models to a sample
predict.consensus.multiclass <-
  function(models,
           data,
           final.class = NA,
           min.prob = 0.25) {
    if (is.null(names(models)))
      names(models) <- 1:length(models)

    lapply(models,
           predict.consensus,
           data = data) %>%
      combine.predicted.list(min.prob = min.prob) ->
      pred

    if (is.na(final.class)) {
      final.class <- paste(names(models), collapse = ".")
    }

    pred %>%
      dplyr::rename_with( ~ paste(.x, final.class, sep = "."), matches("prob$")) %>%
      dplyr::rename_with( ~ final.class, matches("predicted$"))
  }


#' sub-classify a subset of samples belonging to a class
#'
#' returns a table where each rows is predicted outcomes of a sample, by
#' applying a predictive model
#'
#' @param predicted a table where each rows is predicted outcomes of a sample,
#'   i.e. outcome of a \code{predict...} function
#' @param target_column which column to look for class label
#' @param subclassify a class to subclassify
#' @param submodel a list of models, which is an outcome of
#'   \code{\link{cv_loop_train}}
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @returns amended predicted table with subclass labels
predict.consensus.subclass <-
  function(predicted,
           target_column,
           subclassify,
           submodel,
           data) {
    subdata <-
      data[predicted[["Sample_ID"]][predicted[[target_column]] == subclassify],]
    subpredicted <- predict.consensus(submodel, subdata)
    left_join(
      predicted,
      predict.consensus(submodel, subdata) %>%
        dplyr::select("Sample_ID", "prob", "predicted") %>%
        dplyr::rename_with(
          ~ paste("prob", subclassify, "subtype", sep = "."),
          matches("prob$")
        ) %>%
        dplyr::rename_with(
          ~ paste(subclassify, "subtype", sep = "."),
          matches("predicted$")
        ),
      by = "Sample_ID"
    )
  }


#' apply a list of models to a set of samples
#'
#' returns a table where each rows is predicted outcomes of a sample, by
#' applying predictive models
#'
#' @param models a list of models each which is an outcome of
#'   \code{\link{cv_loop_train_iter}} or \code{\link{cv_final_train}}
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param final.class  the name (column heading) of final class outcomes,
#'   majority among multiple classes
#' @param min.prob minimum probability for a sample to be classified to a class
#'   (default: 0.25)
#' @returns a table where each rows is predicted outcomes of applying predictive
#'   models to a sample
predict.type.multiclass <-
  function(models,
           data,
           final.class = NA,
           min.prob = 0.25) {
    if (is.null(names(models)))
      names(models) <- 1:length(models)

    lapply(models,
           predict.type,
           data = data,
           type = "prob") %>%
      combine.predicted.list(min.prob = min.prob) ->
      pred

    if (is.na(final.class)) {
      final.class <- paste(names(models), collapse = ".")
    }

    pred %>%
      dplyr::rename_with( ~ paste(.x, final.class, sep = "."), matches("prob$")) %>%
      dplyr::rename_with( ~ final.class, matches("predicted$"))
  }


#' sub-classify a subset of samples belonging to a class
#'
#' returns a table where each rows is predicted outcomes of a sample, by
#' applying a predictive model
#'
#' @param predicted a table where each rows is predicted outcomes of a sample,
#'   i.e. outcome of a \code{predict...} function
#' @param target_column which column to look for class label
#' @param subclassify a class to subclassify
#' @param submodel a model which is an outcome of
#'   \code{\link{cv_loop_train_iter}} or \code{\link{cv_final_train}}
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @returns amended predicted table with subclass labels
predict.type.subclass <-
  function(predicted,
           target_column,
           subclassify,
           submodel,
           data) {
    subdata <-
      data[predicted[["Sample_ID"]][predicted[[target_column]] == subclassify],]
    subpredicted <- predict.type(submodel, subdata, "prob")
    left_join(
      predicted,
      predict.type(submodel, subdata, "prob") %>%
        dplyr::select("Sample_ID", "prob", "predicted") %>%
        dplyr::rename_with(
          ~ paste("prob", subclassify, "subtype", sep = "."),
          matches("prob$")
        ) %>%
        dplyr::rename_with(
          ~ paste(subclassify, "subtype", sep = "."),
          matches("predicted$")
        ),
      by = "Sample_ID"
    )
  }


#' predict outcomes of observations in data using *selected features*
#'
#' @param model a list of models which is an outcome of
#'   \code{\link{cv_loop_train_iter}} or \code{\link{cv_train_final}}, a model
#'   trained for a class see \code{example} below.
#' @param dd   input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param type 'raw' (default) or 'prob'
#' @return a named vector, either outputs class labels ('raw) or probability
#'   tables, a column for each class label ('prob')
#'
#'          mm = list(
#'            main = a_main_model,
#'            c1 = a_c1_model,
#'            c2 = a_c2_model
#'          )
#'
#'          main is applied first which should classify samples into either c1 or c2
#'          all samples classified to c1 will be then fed into a_c1_model, and
#'          the samples classified to c2 will be fed into a_c2_model.
#'
#'          a model can be NULL, meaning no more action.
#'
#'          The list can be nested.  For example,
#'
#'          mm = list(
#'            main = a_main_model,
#'            c1 = a_c1_model,
#'            c2 = list (
#'              main = a_c2_model,
#'              c2a = a_c2a_model,
#'              c2b = a_c2b_model
#'            )
#'          )
#'
predict.multilayer <- function(models, dd) {
  has.fit <- function(aList) {
    "fit" %in% names(aList)
  }

  classes <- names(models)

  predicted <-
    tibble(Sample_ID = rownames(dd),
           predicted.main = predict.type(models[["main"]], dd, type = "prob"))

  samples_predicted_as_2 <-
    dd[filter(predicted, predicted.main == classes[2])$Sample_ID,]

  samples_predicted_as_3 <-
    dd[filter(predicted, predicted.main == classes[3])$Sample_ID,]

  colnames(predicted)[2] <- paste("pred", classes[1], sep = ".")

  # assume only binary classifiers, so models actually have only three models
  if (has.fit(models[[2]])) {
    predicted_2 <-
      tibble(
        Sample_ID = rownames(samples_predicted_as_2),
        predicted.2 = predict.type(models[[2]], samples_predicted_as_2, type)
      )
    colnames(predicted_2)[2] <- paste("pred", classes[2], sep = ".")

    predicted <- left_join(predicted, predicted_2, by = "Sample_ID")
  } else if (is.list(models[[2]])) {
    predicted_2 <-
      predict.multilayer(models[[2]], samples_predicted_as_2, type)
    colnames(predicted_2)[2] <- paste("pred", classes[2], sep = ".")

    predicted <- left_join(predicted, predicted_2, by = "Sample_ID")
  }

  if (has.fit(models[[3]])) {
    predicted_3 <-
      tibble(
        Sample_ID = rownames(samples_predicted_as_3),
        predicted.3 = predict.type(models[[3]], samples_predicted_as_3, type)
      )
    colnames(predicted_3)[2] <- paste("pred", classes[3], sep = ".")

    predicted <- left_join(predicted, predicted_3, by = "Sample_ID")
  } else if (is.list(models[[3]])) {
    predicted_3 <-
      predict.multilayer(models[[3]], samples_predicted_as_3, type)
    colnames(predicted_3)[2] <- paste("pred", classes[3], sep = ".")

    predicted <- left_join(predicted, predicted_3, by = "Sample_ID")
  }

  predicted
}



#
# PROBABLY OBSOLETE and NO use
#
#' flatten the outcomes of \code{\link{multilayer.predict.type}}
#'
#' @param predicted an outcome of \code{\link{multilayer.predict.type}}
#' @return the last column is added with the multi-class labels (flatted)
flatten.predict.multilayer <- function(predicted) {
  pred <- as.character(predicted[[length(predicted)]])

  for (cidx in seq(length(predicted), 3)) {
    ridx <- which(is.na(pred))
    pred[ridx] <- as.character(predicted[[cidx - 1]][ridx])
  }

  predicted[["predicted"]] <- pred

  predicted
}


#' flatten the outcomes of \code{\link{multilayer.predict.type}}
#'
#' @param predicted an outcome of a \code{.predict....} function
#' @param target_classes a set of classes to be flattened
#' @param final.class  the name (column heading) of final class outcomes, after
#'   flattening classes
#' @return the last column is added with the multi-class labels (flatted)
flatten.predict.multiclass <-
  function(predicted, target_classes, final.class = "final") {
    tc <- predicted[target_classes]
    prob.target_classes <- paste("prob", target_classes, sep = ".")
    i.tc.prob <- match(prob.target_classes, colnames(predicted))
    tc.prob <-
      data.frame(matrix(
        nrow = nrow(predicted),
        ncol = length(target_classes)
      ))
    colnames(tc.prob) <- prob.target_classes
    tc.prob[!is.na(i.tc.prob)] <-
      predicted[prob.target_classes[!is.na(i.tc.prob)]]

    ret <- data.frame(Sample_ID = predicted[["Sample_ID"]],
                      tc.prob, tc)

    pred <- as.character(tc[[length(target_classes)]])
    prob <- tc.prob[[length(prob.target_classes)]]

    for (cidx in seq(length(target_classes), 2)) {
      ridx <- which(is.na(pred))
      pred[ridx] <- as.character(tc[[cidx - 1]][ridx])
      prob[ridx] <- tc.prob[[cidx - 1]][ridx]
    }
    ret[[paste("prob", final.class, sep = ".")]] <- prob
    ret[[final.class]] <- pred
    ret
  }


#' combine outcomes of multiple predictors
#'
#' @param predicted.list a list of predicted outcomes (tables).
#'                       Each table should have two columns: prob, predicted
#' @param min.prob a minimum probability/likelihood to be called a class.
#'                 (default = 0.25)
#' @return combined table
#'         (prob, predicted) of each table will be prefixed with class name,
#'            and combined.  The last two columns will be the final outcomes.
combine.predicted.list <-
  function(predicted.list, min.prob = 0.25) {
    sid <- colnames(predicted.list[[1]])[1]

    pred <- predicted.list[[1]][c(sid, "predicted")]
    for (ii in seq(2, length(predicted.list))) {
      pred <-
        full_join(pred, predicted.list[[ii]][c(sid, "predicted")], by = sid)
    }
    names(pred)[-1] <- names(predicted.list)

    prob <- predicted.list[[1]][c(sid, "prob")]
    for (ii in seq(2, length(predicted.list))) {
      prob <-
        full_join(prob, predicted.list[[ii]][c(sid, "prob")], by = sid)
    }
    names(prob)[-1] <- paste("prob", names(predicted.list), sep = ".")

    combined <-
      data.frame(
        left_join(prob, pred, by = "Sample_ID"),
        idx = 0,
        prob = 0,
        predicted = ""
      )
    for (ii in seq(1, nrow(prob))) {
      idx <- which.max(prob[ii,-1])
      pp <- prob[ii, 1 + idx]
      cls <- colnames(pred)[1 + idx]

      if (pp < min.prob) {
        idx <- which.min(prob[ii,-1])
        cls <- as.character(pred[ii, idx + 1])
      }

      combined[["idx"]][ii] <- idx
      combined[["prob"]][ii] <- pp
      combined[["predicted"]][ii] <- cls
    }

    # list(pred = pred, prob = prob)
    combined
  }




#' predict outcomes of observations in data using *selected features*
#'
#' @param data input matrix, of dimension \code{nobs x nvars};
#'             each row is an observation vector.
#'             Since this is an input to \code{\link{glmnet}},
#'               it should be the format that can be used
#'               with \code{\link{glmnet}}
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
  ret[["trained"]] <- rownames(data) %in% model$train_samples
  ret[["row.names"]] <- rownames(data)
  ret[["pred"]] <- predict.type(model, data, type = "prob")
  # ret[["pred"]] <-
  #   data.frame(pred = predict.type(model, data, type = "raw"),
  #              prob = predict.type(model, data, type = "prob"))
  ret
}



#' Generate a list of prediction summary
#'
#' @param pred_tbls An output from \code{\link{prediction_tables}}
#' @return A list of prediction summary
prediction_summary <- function(pred_tbls) {
  ret <-
    list(
      class = pred_tbls$class,
      trained = pred_tbls$trained,
      row.names = pred_tbls$row.names,
      summary = list()
    )

  ret$summary$p.table <- pred_tbls$pred
  ret$summary$confusion <-
    confusionMatrix(pred_tbls$pred$predicted, pred_tbls$class)
  ret$summary$roc <- pROC::roc(
    pred_tbls$class,
    pred_tbls$pred$prob,
    quiet = TRUE,
    ci = TRUE,
    auc = TRUE,
    plot = FALSE
  )  # prob.[type]

  ret$summary$accuracy <-
    ret$summary$confusion$overall[["Accuracy"]]
  ret$summary$kappa <- ret$summary$confusion$overall[["Kappa"]]
  ret$summary$auc <- ret$summary$roc[["auc"]]
  ret$summary$F1 <- ret$summary$confusion$byClass[["F1"]]
  ret$summary$NPV <-
    ret$summary$confusion$byClass[["Neg Pred Value"]]
  ret$summary$PPV <-
    ret$summary$confusion$byClass[["Pos Pred Value"]]
  ret$summary$sensitivity <-
    ret$summary$confusion$byClass[["Sensitivity"]]
  ret$summary$specificity <-
    ret$summary$confusion$byClass[["Specificity"]]

  ret
}



#' Generate a summary of performance of CV trained models
#'
#' @param cv_trained An output of \code{\link{cv_loop_train}} or
#'   \code{\link{cv_loop_train_parallel}}
models_performance_summary <-
  function(cv_trained) {
    data <- cv_trained$data
    cls <- cv_trained$class
    models <- cv_trained$models

    K <- length(models)

    train_ <- list()
    test_ <- list()
    for (k in 1:K) {
      data_train <- data[models[[k]]$train_index,]
      data_test <- data[-(models[[k]]$train_index),]

      class_train <- cls[models[[k]]$train_index]
      class_test <- cls[-(models[[k]]$train_index)]

      prediction_tables(data = data_train,
                        cls = class_train,
                        model = models[[k]]) %>%
        prediction_summary() -> train_[[k]]

      prediction_tables(data = data_test,
                        cls = class_test,
                        model = models[[k]]) %>%
        prediction_summary() -> test_[[k]]
    }

    list(train = train_, test = test_)
  }



#' Extract performance metrics from \code{performance_summary} slot
#'
#' @param p.summary Either "train" or "test" slot of an output of
#'   \code{\link{models_performance_summary}}
#' @param slot An index
extract_metrics <- function(p.summary) {
  data.frame(
    idx = 1:length(p.summary),
    accuracy = sapply(p.summary, function(x)
      x$summary[["accuracy"]]),
    kappa = sapply(p.summary, function(x)
      x$summary[["kappa"]]),
    # roc = lapply(p.summary, function(x) x$summary[["roc"]]),
    auc = sapply(p.summary, function(x)
      x$summary[["auc"]]),
    F1 = sapply(p.summary, function(x)
      x$summary[["F1"]]),
    NPV = sapply(p.summary, function(x)
      x$summary[["NPV"]]),
    PPV = sapply(p.summary, function(x)
      x$summary[["PPV"]]),
    sensitivity = sapply(p.summary, function(x)
      x$summary[["sensitivity"]]),
    specificity = sapply(p.summary, function(x)
      x$summary[["specificity"]])
  )
}



#
# probably not needed any more
#
#' Extract all performance metrics
#'
#' @param p.summary Either "train" or "test" slot of an output of
#'   \code{\link{models_performance_summary}}
# extract_metrics_all <- function(p.summary) {
#   # p.summary is a array (K trials) of performance summary.
#   # nn gets if both "uniform" and "weighted" are tried.
#   nn <- names(p.summary[[1]]$summary)
#   r <- list()
#   for (n in nn) {
#     r[[n]] <- extract_metrics(p.summary, slot = n)
#   }
#
#   bind_rows(r, .id = "feature_weight")
# }


#' Construct feature maps
#'
#' @param models a list of models each of which is an output of
#'   \code{\link{cv_loop_train_iter}}
#' @return a list ("uniform" and/or "weighted") of feature map which is made of
#'   * columns of filtered features, amended by \code{overall_score} *
#'   predictors (selected by model), amended by \code{overall_score,
#'   overall_feature_score} * model coefficients, amended by
#'   \code{overall_score, overall_score_mad, overall_score2,
#'   overall_feature_score}
construct_feature_maps <- function(models) {
  # features (filtered)
  lapply(models, function(m)
    m$selected_features) %>%
    bind_rows(.id = "id") %>%
    dplyr::select(id, feature, score) %>%
    pivot_wider(
      names_from = "id",
      values_from = "score",
      values_fill = 0
    ) %>%
    data.frame -> feature_map

  # predictors (selected by model)
  lapply(models, function(m)
    m$selected_features) %>%
    bind_rows(.id = "id") %>%
    dplyr::select(id, feature, predictor) %>%
    pivot_wider(
      names_from = "id",
      values_from = "predictor",
      values_fill = 0
    ) %>%
    data.frame -> predictor_map

  # model coefficients
  lapply(models, function(m) {
    coef(m$fit$finalModel,
         s = m$fit$bestTune$lambda) %>% as.matrix %>% as.data.frame -> x
    colnames(x)[1] <- "coef"
    x[["feature"]] <- rownames(x)
    rownames(x) <- NULL
    x
  }) %>%
    bind_rows(.id = "id") %>%
    #      filter(feature != "(Intercept)") %>%
    dplyr::select(id, feature, coef) %>%
    pivot_wider(
      names_from = "id",
      values_from = "coef",
      values_fill = 0
    ) %>%
    data.frame -> coef_map

  rownames(feature_map) <- feature_map$feature
  rownames(predictor_map) <- predictor_map$feature

  # set aside "(Intercept)"
  coef_intercept <- filter(coef_map, feature == "(Intercept)")
  coef_map <- filter(coef_map, feature != "(Intercept)")

  # remove those features never selected
  rx <- rowSums(predictor_map[-1])
  predictor_map <- predictor_map[rx > 0,]

  rx <- rowSums(abs(coef_map[-1]))
  coef_map <- coef_map[rx > 0,]

  # feature scores (overall), from individual scores
  feature_map[["overall_score"]] <-
    feature_map %>%
    dplyr::select(starts_with("X")) %>%
    rowMeans()

  # predictor score = # of times selected in a model
  predictor_map[["overall_score"]] <-
    predictor_map %>%
    dplyr::select(starts_with("X")) %>%
    rowSums()

  f2p_mapping <- match(predictor_map$feature, feature_map$feature)
  predictor_map[["overall_feature_score"]] <-
    feature_map[["overall_score"]][f2p_mapping]

  # overall contribution positive or negative.
  coef_map[["overall_score"]] <-
    coef_map %>%
    dplyr::select(starts_with("X")) %>%
    sign () %>% rowMeans()

  coef_map[["overall_score_mad"]] <-
    coef_map %>%
    dplyr::select(starts_with("X")) %>%
    as.matrix() %>% rowMads()

  coef_map[["overall_score2"]] <-
    coef_map[["overall_score_mad"]] / coef_map[["overall_score"]]

  f2p_mapping <- match(coef_map$feature, feature_map$feature)

  coef_map[["overall_feature_score"]] <-
    feature_map[["overall_score"]][f2p_mapping]


  feature_map %>%
    arrange(-overall_score) %>%
    dplyr::select(feature, overall_score, starts_with("X")) ->
    feature_map

  predictor_map %>%
    arrange(-overall_score) %>%
    dplyr::select(feature,
                  overall_score,
                  overall_feature_score,
                  starts_with("X")) ->
    predictor_map

  coef_map %>%
    arrange(-overall_score) %>%
    dplyr::select(feature,
                  starts_with("overall_score"),
                  overall_feature_score,
                  starts_with("X")) ->
    coef_map

  list(
    feature_map = feature_map,
    predictor_map = predictor_map,
    coef_map = coef_map,
    coef_intercept = coef_intercept
  )
}


#' Workflow to get all classification performance
#'
#' @param cv_trained An output of \code{\link{cv_loop_train}} or
#'   \code{\link{cv_loop_train_parallel}}
#' @return A list of classification performance summary
classification_summary_workflow <-
  function(cv_trained, show_feature_map = FALSE) {
    cv_trained_summary <- list(cv_trained = cv_trained)

    cv_trained_summary[["performance_summary"]] <-
      models_performance_summary(cv_trained)
    cv_trained_summary[["accuracy"]] <-
      bind_rows(
        cbind(
          feature_weight = cv_trained$models[[1]]$feature_weights,
          type = "train",
          extract_metrics(cv_trained_summary[["performance_summary"]][["train"]])
        ),
        cbind(
          feature_weight = cv_trained$models[[1]]$feature_weights,
          type = "test",
          extract_metrics(cv_trained_summary[["performance_summary"]][["test"]])
        )
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
    cv_trained_summary[["accuracy_violin"]] <-
      wrap_up_accuracy_gp (accuracy_gp +
                             geom_violin(aes(color = name)))


    # box plots
    cv_trained_summary[["accuracy_boxplot"]] <-
      wrap_up_accuracy_gp (accuracy_gp +
                             geom_boxplot(aes(color = name), notch = TRUE, width = 0.8))


    # print(cv_trained_summary[["accuracy_boxplot"]])

    # summary features
    cv_trained_summary[["feature_maps"]] <-
      construct_feature_maps(cv_trained$models)

    cv_trained_summary <-
      construct_feature_heatmaps(cv_trained_summary)

    # temp variable to ease the coding below
    feature_maps <- cv_trained_summary[["feature_maps"]]

    # count the number of features/predictors
    n_features <- function(x_map) {
      colSums(dplyr::select(x_map, starts_with("X")) > 0)
    }

    n_features_list <- list()
    if (!is.null(feature_maps[["uniform"]])) {
      n_features_list[["uniform.n_features"]] <-
        n_features(feature_maps[["uniform"]]$feature_map)
      n_features_list[["uniform.n_predictors"]] <-
        n_features(feature_maps[["uniform"]]$predictor_map)

      if (show_feature_map)
        print(feature_maps[["uniform"]]$coef_heatmap)

    }

    if (!is.null(feature_maps[["weighted"]])) {
      n_features_list[["weighted.n_features"]] <-
        n_features(feature_maps[["weighted"]]$feature_map)
      n_features_list[["weighted.n_predictors"]] <-
        n_features(feature_maps[["weighted"]]$predictor_map)

      if (show_feature_map)
        print(feature_maps[["weighted"]]$coef_heatmap)
    }

    cv_trained_summary[["n_features"]] <- n_features_list

    cv_trained_summary
  }



#' Generate a consensus summary of outcomes
#'
#' @param cv_trained An output of \code{\link{cv_loop_train}} or
#'   \code{\link{cv_loop_train_parallel}}
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
    x2[[cls]] <- rowSums(x2[,-c(1:2)] == cls)
  }

  x2[["majority"]] <-
    factor(cls_labels[apply(x2[cls_labels], MARGIN = 1, which.max)], cls_labels)
  x2[["purity"]] <-
    apply(x2[cls_labels], MARGIN = 1, max) / length(cv_trained$models)

  x2 <- mutate(x2, mistake = cls != majority)

  x2 %>%
    #    select(c(1:2), cls_labels, majority, purity, mistake, starts_with("X")) %>%
    dplyr::select(-starts_with("X"), starts_with("X")) %>%
    arrange(-mistake, purity)
}
