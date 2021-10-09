#' A function to compute difference between two classes
#'
#' @param dd A \code{matrix} or \code{data.frame}
#'           where each row is an observation
#' @param cls A set of classes, either in \code{factor}
#' @return A vector of difference
#'            \code{dd[cls == factor.1] - dd[cls==factor.2]}
compute_difference <- function(dd, cls) {
  # to assign the score with direction
  lvls <- levels(cls)
  if (length(lvls) > 2) {
    stop("class should be binary...")
  }

  colMeans(dd[cls == lvls[1],]) -
    colMeans(dd[cls == lvls[2],])
}



#' A function for filtering-based feature selection, using ROC/AUC as a metric
#'
#' @param dd A \code{matrix} or \code{data.frame}
#'           where each row is an observation
#' @param cls A set of classes, either in \code{factor}
#' @param threshold A minimum value of AUC for filtering
#' @param direction Direction of AUC statistical test: "two.sided", "greater" or "less"
#' @return A \code{data.frame} of \code{feature}, \code{diff}, \code{score}
#'         where \code{score} > \code{threshold}
filter_features_AUC <-
  function(dd,
           cls,
           threshold = 0.8,
           direction = c("two.sided", "greater", "less")) {
    res <- caret::filterVarImp(x = dd, y = cls)
    res <-
      dplyr::bind_cols(tibble::tibble(feature = rownames(res)),
                       diff = compute_difference(dd, cls),
                       score = res[[1]])

    dir <- pmatch(direction[1], c("two.sided", "greater", "less"))

    if (dir == 1) {
      res <-
        dplyr::arrange(res,-score) %>% dplyr::filter(score > threshold)
    } else if (dir == 2) {
      res <-
        dplyr::arrange(res,-score) %>% dplyr::filter(diff > 0, score > threshold)
    } else if (dir == 3) {
      res <-
        dplyr::arrange(res,-score) %>% dplyr::filter(diff < 0, score > threshold)
    }

    res
  }

#' Internal function for multiple Wilcox test
#'
#'   only works for binary
#'   cls should be 0 or 1
#'
#' @param dd A \code{matrix} or \code{data.frame}
#'           where each row is an observation
#' @param cls A set of classes, either in \code{factor}
mult_wilcox_test <- function(dd, cls, ...) {
  x <- apply(dd, MARGIN = 1,
             function(x)
               stats::wilcox.test(x[cls == 0], x[cls == 1], ...))
  data.frame(
    p.value = sapply(x, function(xt)
      xt$p.value),
    stats = sapply(x, function(xt)
      xt$statistic)
  )
}



#' A function for filtering-based feature selection, using Wilcox test as a metric
#'
#' @param dd A \code{matrix} or \code{data.frame}
#'           where each row is an observation
#' @param cls A set of classes, either in \code{factor}
#' @param threshold A minimum value of (1 - p.value of Wilcox test) for filtering
#' @param direction Direction of AUC statistical test: "two.sided", "greater" or "less"
#' @return A table where each row is a filtered feature with stats and scores
filter_features <-
  function(dd,
           cls,
           threshold = 0.8,
           direction = c("two.sided", "greater", "less")) {
    res <-
      mult_wilcox_test(t(dd), as.numeric(cls) - min(as.numeric(cls)), alternative = direction[1])
    res <- dplyr::bind_cols(
      tibble::tibble(feature = rownames(res)),
      res,
      diff = compute_difference(dd, cls),
      score = 1 - res$p.value
    )

    dir <- pmatch(direction[1], c("two.sided", "greater", "less"))

    if (dir == 1) {
      res <-
        dplyr::arrange(res,-score) %>% dplyr::filter(score > threshold)
    } else if (dir == 2) {
      res <-
        dplyr::arrange(res,-score) %>% dplyr::filter(diff > 0, score > threshold)
    } else if (dir == 3) {
      res <-
        dplyr::arrange(res,-score) %>% dplyr::filter(diff < 0, score > threshold)
    }

    res
  }



#' A function train a model with CV
#'
#' @param data_train input matrix, of dimension nobs x nvars; each row is an
#'   observation vector. Since this is an input to \code{glmnet}, it should be
#'   the format that can be used with \code{glmnet}
#' @param cls_train class labels
#' @param fitControl A list of training parameters.  See
#'   \code{\link{caret::trainControl}} for detail
#' @param penalty.factor Separate penalty factors can be applied to each
#'   coefficient.  See \code{\link{glmnet::glmnet}} for detail
#' @return A list returned from \code{caret::train}
train_cv_model_glmnet <- function(data_train,
                                  cls_train,
                                  fitControl,
                                  observation_weights = NULL,
                                  penalty.factor = 1) {
  if (length(penalty.factor) == 1)
    penalty.factor <- rep(1, ncol(data_train))

  # Since `lambda` is critical in glmnet (lasso, alpha = 1), let's get initial grid for that.
  # We will do first this with penalty.factor
  #
  glmnet_fit <- glmnet::glmnet(
    x = as.matrix(data_train),
    y = cls_train,
    weights = observation_weights,
    penalty.factor = penalty.factor,
    family = "binomial"
  )

  # Only Lasso, hence, alpha = 1
  tuneGrid = expand.grid(.alpha = 1, .lambda = glmnet_fit$lambda)

  #
  # actual training with CV built-in, which will be done according to 'fitControl'
  #   fitControl defines # of CV folds and # of CV repeats
  #
  return(
    caret::train(
      x = data_train,
      y = cls_train,
      method = "glmnet",
      weights = observation_weights,
      penalty.factor = penalty.factor,
      trControl = fitControl,
      tuneGrid = tuneGrid
    )
  )
}


#' A function to get selected features from a model, \code{fit}
#'
#' @param fit ...
#' @param features ...
#' @return A table where each feature (row) is marked 1 if it is a selected
#'   feature, 0 otherwise.
get_selected_features <- function(fit, features) {
  selected_predictors <-
    tibble::tibble(feature = caret::predictors(fit),
                   predictor = 1)

  dplyr::left_join(features,
                   selected_predictors, by = "feature") %>%
    dplyr::mutate(predictor = ifelse(is.na(predictor), 0, predictor))
}


#' A function to train a model with CV
#'
#' @param ii integer (probably not needed)
#' @param data input matrix, of dimension nobs x nvars; each row is an
#'   observation vector. Since this is an input to \code{glmnet}, it should be
#'   the format that can be used with \code{glmnet}
#' @param cls class labels
#' @param fitControl A list of training parameters.  See
#'   \code{\link{caret::trainControl}} for detail
#' @param resampling_rate ...
#' @param n_features ...
#' @param filter_method ...
#' @param filter_threshold_score ...
#' @param filter_threshold_diff ...
#' @param filter_direction ...
#' @param feature_uniform ...
#' @param feature_weighted ...
#' @return a trained model - a list of training parameters (see above), fitted
#'   model, and selected features
cv_loop_train_iter <-
  function(ii,
           data,
           cls,
           fitControl,
           resampling_rate,
           n_features,
           filter_method,
           filter_threshold_score,
           filter_threshold_diff,
           filter_direction,
           observation_weights,
           feature_weights) {
    aModel <- list()

    stopifnot(is.factor(cls))

    maxTried = 5  # it was 20 initially, but reduced for computation and
    #     because it won't make sense anyway if you don't get
    #     enough filtered features within a few trials.
    for (nTried in 0:maxTried) {
      train_index <-
        caret::createDataPartition(cls,
                                   p = resampling_rate,
                                   list = FALSE,
                                   times = 1)

      data_train <- data[train_index,]
      cls_train <- cls[train_index]

      if (!is.null(observation_weights)) {
        observation_weights_train <- observation_weights[train_index]
      } else {
        observation_weights_train <- observation_weights
      }

      # Filtering-based feature selection, candidates for the wrapper-based feature selection
      if (toupper(filter_method) == "WILCOX") {
        filtered_features <-
          filter_features(
            dd = data_train,
            cls = cls_train,
            threshold = 0,
            direction = filter_direction
          )
      } else if (toupper(filter_method) == "ROC") {
        filtered_features <-
          filter_features_AUC(
            dd = data_train,
            cls = cls_train,
            threshold = 0,
            direction = filter_direction
          )
      } else {
        stop("'filter_method' should be either 'ROC' or 'wilcox'.")
      }

      filtered_features <-
        dplyr::filter(filtered_features, abs(diff) > filter_threshold_diff)

      #
      # if n_features is not pre-defined,
      #   we will need filter_threshold_score to pick a subset of features for the next step
      #
      if (is.na(n_features)) {
        filtered_features <-
          dplyr::filter(filtered_features, score > filter_threshold_score)
      } else {
        filtered_features %>%
          dplyr::arrange(-score) %>%
          dplyr::slice_head(n = n_features) ->
          filtered_features
      }

      # somehow, the number of selected_features is less than 2, we need to repeat the step above.
      if (nrow(filtered_features) > 1)
        break
    }

    if (nrow(filtered_features) < 2) {
      stop("# of filtered_features is less than 2; try to loosen the filtering criteria...")
    }

    # we will need data with only selected features
    data_train_with_filtered_features <-
      data_train[, filtered_features$feature]

    aModel <- list(
      train_index = train_index,
      train_samples = rownames(data_train_with_filtered_features),
      filter_method = filter_method[1],
      filter_threshold_diff = filter_threshold_diff,
      filter_threshold_score = filter_threshold_score,
      filter_direction = filter_direction,
      observation_weights = observation_weights_train,
      feature_weights = feature_weights
    )

    if (feature_weights == "weighted") {
      penalty.factor <- (1 - filtered_features$score)
    } else {
      penalty.factor <- 1
    }


    fit_cv <-
      train_cv_model_glmnet(
        data_train_with_filtered_features,
        cls_train,
        fitControl,
        observation_weights = observation_weights_train,
        penalty.factor = penalty.factor
      )
    selected_features <-
      get_selected_features(fit_cv, filtered_features)

    aModel[["fit"]] <- fit_cv
    aModel[["selected_features"]] <- selected_features

    aModel
  }




#' iterates training a model with CV (serial version)
#'
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param cls class labels
#' @param K ...
#' @param resampling_rate ...
#' @param n_features ...
#' @param filter_method ...
#' @param filter_threshold_score ...
#' @param filter_threshold_diff ...
#' @param filter_direction ...
#' @param observation_weights ...
#' @param feature_weights ...
#' @return a list of trained models -- see \code{\link{cv_loop_train_iter}}
cv_loop_train <- function(data,
                          cls,
                          fitControl,
                          K = 25,
                          resampling_rate = 0.8,
                          n_features = NA,
                          filter_method = c("ROC", "WILCOX"),
                          filter_direction = c("two.sided", "greater", "less"),
                          filter_threshold_diff = 1,
                          filter_threshold_score = 0.8,
                          observation_weights = NULL,
                          feature_weights = c("uniform", "weighted")) {
  k <- pmatch(feature_weights, c("uniform", "weighted"))
  if (is.na(k)) {
    stop("feature_weights should be either 'uniform' or 'weighted'.")
  }

  feature_weights <- feature_weights[k]

  cv_model_list <-
    purrr::map(
      1:K,
      cv_loop_train_iter,
      data = data,
      cls = cls,
      fitControl = fitControl,
      resampling_rate = resampling_rate,
      n_features = n_features,
      filter_method = filter_method[1],
      filter_direction = filter_direction[1],
      filter_threshold_diff = filter_threshold_diff,
      filter_threshold_score = filter_threshold_score,
      observation_weights = observation_weights,
      feature_weights = feature_weights[1]
    )

  list(data = data,
       class = cls,
       models = cv_model_list)
}


#' iterates training a model with CV (parallel version), using
#' \code{\link{future.apply}}
#'
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an
#'   observation vector. Since this is an input to \code{\link{glmnet}}, it
#'   should be the format that can be used with \code{\link{glmnet}}
#' @param cls class labels
#' @param K ...
#' @param resampling_rate ...
#' @param n_features ...
#' @param filter_method ...
#' @param filter_threshold_score ...
#' @param filter_threshold_diff ...
#' @param filter_direction ...
#' @param observation_weights ...
#' @param feature_weights ...
#' @return a list of trained models -- see
#'   \code{\link{CCSBUtil::cv_loop_train_iter}}
cv_loop_train_parallel <-
  function(data,
           cls,
           fitControl,
           K = 25,
           resampling_rate = 0.8,
           n_features = NA,
           filter_method = c("ROC", "WILCOX"),
           filter_direction = c("two.sided", "greater", "less"),
           filter_threshold_diff = 1,
           filter_threshold_score = 0.8,
           observation_weights = NULL,
           feature_weights = c("uniform", "weighted")) {
    k <- pmatch(feature_weights, c("uniform", "weighted"))
    if (is.na(k)) {
      stop("feature_weights should be either 'uniform' or 'weighted'.")
    }

    feature_weights <- feature_weights[k]

    cv_model_list <-
      future.apply::future_lapply(
        1:K,
        cv_loop_train_iter,
        data = data,
        cls = cls,
        fitControl = fitControl,
        resampling_rate = resampling_rate,
        n_features = n_features,
        filter_method = filter_method[1],
        filter_direction = filter_direction[1],
        filter_threshold_diff = filter_threshold_diff,
        filter_threshold_score = filter_threshold_score,
        observation_weights = observation_weights,
        feature_weights = feature_weights[1],
        future.seed = TRUE
      )

    list(data = data,
         class = cls,
         models = cv_model_list)
  }


#' A function to train a *final* model with CV
#'
#' @param data input matrix, of dimension nobs x nvars; each row is an
#'   observation vector. Since this is an input to \code{glmnet}, it should be
#'   the format that can be used with \code{glmnet}
#' @param cls class labels
#' @param fitControl A list of training parameters.  See
#'   \code{\link{caret::trainControl}} for detail
#' @param filter_method ... c("ROC", "WILCOX")
#' @param filter_direction ... c("two.sided", "greater", "less")
#' @param observation_weights ...
#' @param feature_weights ... c("uniform", "weighted")
#' @return a trained model (final) - a list of training parameters (see above),
#'   fitted model, and selected features
cv_train_final <- function(data,
                           cls,
                           fitControl,
                           filter_method = c("ROC", "WILCOX"),
                           filter_direction = c("two.sided", "greater", "less"),
                           observation_weights = NULL,
                           feature_weights = c("uniform", "weighted")) {
  aModel <- list()

  filter_method <- filter_method[1]
  filter_direction <- filter_direction[1]
  feature_weights <- feature_weights[1]

  # Filtering-based feature selection, candidates for the wrapper-based feature selection
  if (toupper(filter_method) == "WILCOX") {
    filtered_features <-
      filter_features(
        dd = data,
        cls = cls,
        threshold = 0,
        direction = filter_direction
      )
  } else if (toupper(filter_method) == "ROC") {
    filtered_features <-
      filter_features_AUC(
        dd = as.data.frame(data),
        # this should work w/o as.data.frame. should check it out later.
        cls = cls,
        threshold = 0,
        direction = filter_direction
      )
  } else {
    stop("'filter_method' should be either 'ROC' or 'wilcox'.")
  }

  filtered_features <-
    dplyr::filter(filtered_features, abs(diff) > 0)

  aModel <- list(
    train_index = 1:nrow(data),
    train_samples = rownames(data),
    filter_method = filter_method[1],
    filter_threshold_diff = NA,
    filter_threshold_score = NA,
    filter_direction = filter_direction[1],
    observation_weights = observation_weights,
    feature_weights = feature_weights
  )

  if (feature_weights == "weighted") {
    penalty.factor <- (1 - filtered_features$score)
  } else {
    penalty.factor <- 1
  }

  # making sure the features are in the right order
  data <- data[, filtered_features$feature]
  fit_cv <- train_cv_model_glmnet(data,
                                  cls,
                                  fitControl,
                                  observation_weights = observation_weights,
                                  penalty.factor = penalty.factor)

  selected_features <-
    get_selected_features(fit_cv, filtered_features)

  aModel[["fit"]] <- fit_cv
  aModel[["selected_features"]] <- selected_features

  aModel
}


#' A function that combines cv_loop_train and cv_train_final
#'
#' @param data input matrix, of dimension nobs x nvars; each row is an
#'   observation vector. Since this is an input to \code{glmnet}, it should be
#'   the format that can be used with \code{glmnet}
#' @param cls class labels
#' @param fitControl A list of training parameters.  See
#'   \code{\link{caret::trainControl}} for detail
#' @param K ...
#' @param resampling_rate ...
#' @param n_features ...
#' @param filter_method ...
#' @param filter_threshold_score ...
#' @param filter_threshold_diff ...
#' @param filter_direction ...
#' @param observation_weights ...
#' @param feature_weights ...
#' @param predictor_score_threshold ...
#' @return a list of \code{cv_loop_trained} (see
#'   \code{\link{cv_loop_train_iter}}), \code{classification_results} (see
#'   \code{\link{classification_summary_workflow}}), and \code{final_cv_model}
#'   (see \code{\link{cv_train_final}} and additional slots)
train_to_final_model <-
  function(data,
           cls,
           train_final = TRUE,
           train_consensus = TRUE,
           fitControl,
           K = 25,
           resampling_rate = 0.8,
           n_features = NA,
           filter_method = c("ROC", "WILCOX"),
           filter_direction = c("two.sided", "greater", "less"),
           filter_threshold_diff = 1,
           filter_threshold_score = 0.8,
           observation_weights = NULL,
           feature_weights = c("uniform", "weighted"),
           predictor_score_threshold = 0.1,
           verbose = FALSE) {
    if (train_consensus) {
      #
      # cv_loop training to get consensus model
      #
      cv_loop_trained <-
        cv_loop_train_parallel(
          data = data,
          cls = cls,
          fitControl = fitControl,
          K = K,
          resampling_rate = resampling_rate,
          filter_method = filter_method,
          filter_direction = filter_direction,
          filter_threshold_diff = filter_threshold_diff,
          filter_threshold_score = filter_threshold_score,
          n_features = n_features,
          feature_weights = feature_weights
        )

      if (verbose) {
        print("cv_loop_train_parallel complete ...")
      }

      cv_loop_trained %>%
        classification_summary_workflow() ->
        classification_results

      cv_loop_trained %>%
        prediction_consensus_summary() ->
        classification_results[["consensus"]]

      classification_results$feature_maps$coef_map %>%
        filter(abs(overall_score) > predictor_score_threshold) %>%
        `[[`("feature") ->
        features_selected

      features_df <-
        enframe(
          compute_difference(
            classification_results$cv_trained$data[, features_selected],
            classification_results$cv_trained$class
          )
        ) %>%
        dplyr::rename(diff = value)

      classification_results[["feature_boxplots"]] <-
        feature_boxplots(
          dd = classification_results$cv_trained$data[, features_selected],
          cls = classification_results$cv_trained$class,
          df = features_df,
          order.by = "diff"
        )

      if (verbose) {
        print("classification_summary_workflow complete ...")
      }

    } else {
      cv_loop_trained <- NA
      classification_results <- NA
      features_selected <- rownames(data)  # use all features
    }

    if (verbose) {
      print("Consensus model complete ...")
    }

    if (train_final) {
      final_data <- data[features_selected]
      final_cls <- cls

      final_cv_model <-
        cv_train_final(
          data = final_data,
          cls = final_cls,
          fitControl = fitControl,
          filter_method = filter_method
        )

      final_cv_model[["summary"]] <-
        data.frame(
          sample = rownames(final_cv_model$fit$trainingData),
          prob = predict(final_cv_model$fit, type = "prob"),
          pred = predict(final_cv_model$fit),
          cls = final_cv_model$fit$trainingData$.outcome
        ) %>%
        mutate(correct = cls == pred) %>%
        arrange(-correct)

      final_features_selected <- predictors(final_cv_model$fit)

      final_features_df <-
        left_join(enframe(compute_difference(final_data[, final_features_selected],
                                             final_cls)),
                  as_tibble(varImp(final_cv_model$fit)[[1]],
                            rownames = 'name'),
                  by = 'name') %>%
        dplyr::rename(diff = value,
                      imp = Overall)

      final_cv_model[["candidate_features"]] <- features_selected
      final_cv_model[["final_features"]] <- final_features_selected

      final_cv_model[["summary"]] %>%
        dplyr::select(starts_with("prob.")) %>%
        colnames ->
        cls.probs

      final_cv_model[["plots"]] <-
        list(
          final_train_trajectory = ggplot_train_trajectory(final_cv_model$fit),
          final_prediction_probability = final_cv_model[["summary"]] %>%
            ggplot(aes(
              x = reorder(sample, eval(parse(text = cls.probs[1]))),
              y = eval(parse(text = cls.probs[1])),
              color = cls
            )) +
            geom_point() +
            ylab(cls.probs[1]) +
            xlab("samples") +
            geom_hline(
              yintercept = 0.5,
              linetype = "dashed",
              color = "red"
            ) +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()),
          final_features_importance = ggplot(varImp(final_cv_model$fit)),
          final_features_boxplots = feature_boxplots(
            dd = final_data[, final_features_selected],
            cls = final_cls,
            df = final_features_df,
            order.by = "diff"
          )
        )
    } else {
      final_cv_model <- NA
    }

    if (verbose) {
      print("Final model complete ...")
    }

    list(
      cv_loop_trained = cv_loop_trained,
      classification_results = classification_results,
      final_cv_model = final_cv_model
    )
  }
