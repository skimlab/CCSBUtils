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

  colMeans(dd[cls == lvls[1], ]) -
    colMeans(dd[cls == lvls[2], ])
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
           direction = c("two.sided", "greater", "less")
  ) {
  res <- caret::filterVarImp(x = dd, y = cls)
  res <-
    dplyr::bind_cols(tibble::tibble(feature = rownames(res)),
                     diff = compute_difference(dd, cls),
                     score = res[[1]])

  dir <- pmatch(direction, c("two.sided", "greater", "less"))[1]

  if (dir == 1) {
    res <- dplyr::arrange(res, -score) %>% dplyr::filter(score > threshold)
  } else if (dir == 2) {
    res <- dplyr::arrange(res, -score) %>% dplyr::filter(diff > 0, score > threshold)
  } else if (dir == 3) {
    res <- dplyr::arrange(res, -score) %>% dplyr::filter(diff < 0, score > threshold)
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
               stats::wilcox.test(x[cls==0], x[cls==1], ...))
  data.frame(p.value = sapply(x, function(xt) xt$p.value), stats = sapply(x, function(xt) xt$statistic))
}



#' A function for filtering-based feature selection, using Wilcox test as a metric
#'
#' @param dd A \code{matrix} or \code{data.frame}
#'           where each row is an observation
#' @param cls A set of classes, either in \code{factor}
#' @param threshold A minimum value of (1 - p.value of Wilcox test) for filtering
#' @param direction Direction of AUC statistical test: "two.sided", "greater" or "less"
filter_features <-
  function(dd,
           cls,
           threshold = 0.8,
           direction = c("two.sided", "greater", "less")
  ) {
  res <- mult_wilcox_test(t(dd), as.numeric(cls)-1, alternative = direction[1])
  res <- dplyr::bind_cols(tibble::tibble(feature = rownames(res)),
                          res,
                          diff = compute_difference(dd, cls),
                          score = 1 - res$p.value)

  dir <- pmatch(direction, c("two.sided", "greater", "less"))[1]

  if (dir == 1) {
    res <- dplyr::arrange(res, -score) %>% dplyr::filter(score > threshold)
  } else if (dir == 2) {
    res <- dplyr::arrange(res, -score) %>% dplyr::filter(diff > 0, score > threshold)
  } else if (dir == 3) {
    res <- dplyr::arrange(res, -score) %>% dplyr::filter(diff < 0, score > threshold)
  }

  res
}



#' A function train a model with CV
#'
#' @param data_train input matrix, of dimension nobs x nvars; each row is an observation vector.
#'                   Since this is an input to \code{glmnet}, it should be the format that can be used
#'                   with \code{glmnet}
#' @param cls_train class labels
#' @param fitControl A list of training parameters.  See \code{\link{caret::trainControl}} for detail
#' @param penalty.factor Separate penalty factors can be applied to each coefficient.  See \code{\link{glmnet::glmnet}} for detail
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
      tuneGrid = tuneGrid)
  )
}


#' A function to get selected features from a model, \code{fit}
#'
#' @param fit ...
#' @param features ...
#'
get_selected_features <- function(fit, features) {
  selected_predictors <-
    tibble::tibble(feature = caret::predictors(fit),
                   predictor = 1)

  dplyr::left_join(features,
            selected_predictors, by = "feature") %>%
    dplyr::mutate(predictor = ifelse(is.na(predictor), 0, predictor))
}


#' A function train a model with CV
#'
#' @param ii integer (probably not needed)
#' @param data input matrix, of dimension nobs x nvars; each row is an observation vector.
#'             Since this is an input to \code{glmnet}, it should be the format that can be used
#'             with \code{glmnet}
#' @param cls class labels
#' @param fitControl A list of training parameters.  See \code{\link{caret::trainControl}} for detail
#' @param resampling_rate ...
#' @param n_features ...
#' @param filter_method ...
#' @param filter_threshold_score ...
#' @param filter_threshold_diff ...
#' @param filter_direction ...
#' @param feature_uniform ...
#' @param feature_weighted ...
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

    repeat {
      train_index <- caret::createDataPartition(cls, p = resampling_rate, list = FALSE, times = 1)

      data_train <- data[train_index, ]
      cls_train <- cls[train_index]

      if (!is.null(observation_weights)) {
        observation_weights_train <- observation_weights[train_index]
      } else {
        observation_weights_train <- observation_weights
      }

      # Filtering-based feature selection, candidates for the wrapper-based feature selection
      if (toupper(filter_method) == "WILCOX") {
        filtered_features <-
          filter_features(dd = data_train, cls = cls_train,
                          threshold = 0, direction = filter_direction)
      } else if (toupper(filter_method) == "ROC") {
        filtered_features <-
          filter_features_AUC(dd = data_train, cls = cls_train,
                              threshold = 0, direction = filter_direction)
      } else {
        stop("'filter_method' should be either 'ROC' or 'wilcox'.")
      }

      filtered_features <- dplyr::filter(filtered_features, abs(diff) > filter_threshold_diff)

      #
      # if n_features is not pre-defined,
      #   we will need filter_threshold_score to pick a subset of features for the next step
      #
      if (is.na(n_features)) {
        filtered_features <- dplyr::filter(filtered_features, score > filter_threshold_score)
      } else {
        filtered_features %>%
          dplyr::arrange(-score) %>%
          dplyr::slice_head(n = n_features) ->
          filtered_features
      }

      # somehow, the number of selected_features is 0, we need to repeat the step above.
      if (nrow(filtered_features) > 0)
        break
    }

    # we will need data with only selected features
    data_train_with_filtered_features <- dplyr::select(data_train, filtered_features$feature)

    aModel <- list(
      train_index = train_index,
      train_samples = rownames(data_train_with_filtered_features),
      filter_method = filter_method,
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


    fit_cv <- train_cv_model_glmnet(data_train_with_filtered_features,
                                    cls_train,
                                    fitControl,
                                    observation_weights = observation_weights_train,
                                    penalty.factor = penalty.factor)
    selected_features <- get_selected_features(fit_cv, filtered_features)

    aModel[["fit"]] <- fit_cv
    aModel[["selected_features"]] <- selected_features

    aModel
  }


#' iterate training a model with CV
#'
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an observation vector.
#'             Since this is an input to \code{\link{glmnet}}, it should be the format that can be used
#'             with \code{\link{glmnet}}
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
cv_loop_train <- function(data,
                          cls,
                          fitControl,
                          K = 25,
                          resampling_rate = 0.8,
                          n_features = NA,
                          filter_method = c("two.sided", "greater", "less"),
                          filter_direction = "less",
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
      filter_method = filter_method,
      filter_direction = filter_direction,
      filter_threshold_diff = filter_threshold_diff,
      filter_threshold_score = filter_threshold_score,
      observation_weights = observation_weights,
      feature_weights = feature_weights
    )

  list(data = data, class = cls, models = cv_model_list)
}


#' iterate training a model with CV, using \code{\link{future.apply}}
#'
#' @param data input matrix, of dimension \code{nobs x nvars}; each row is an observation vector.
#'             Since this is an input to \code{\link{glmnet}}, it should be the format that can be used
#'             with \code{\link{glmnet}}
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
cv_loop_train_parallel <-
  function(data,
           cls,
           fitControl,
           K = 25,
           resampling_rate = 0.8,
           n_features = NA,
           filter_method = c("two.sided", "greater", "less"),
           filter_direction = "less",
           filter_threshold_diff = 1,
           filter_threshold_score = 0.8,
           observation_weights = NULL,
           feature_weights = c("uniform", "weighted")
           ) {

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
        filter_method = filter_method,
        filter_direction = filter_direction,
        filter_threshold_diff = filter_threshold_diff,
        filter_threshold_score = filter_threshold_score,
        observation_weights = observation_weights,
        feature_weights = feature_weights,
        future.seed = TRUE
      )

    list(data = data, class = cls, models = cv_model_list)
  }



