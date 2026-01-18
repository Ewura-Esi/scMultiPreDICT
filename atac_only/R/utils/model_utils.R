# ============================================================================
# scMultiPreDICT - Model Training Utilities
# ============================================================================
# Reusable functions for training predictive models
#
# Usage:
#   source("R/utils/model_utils.R")
# ============================================================================

#' Train all linear models for a target gene
#'
#' @param X Feature matrix (cells × features)
#' @param y Target expression vector
#' @param alpha_values Alpha values for elastic net (0=ridge, 1=lasso)
#' @param lambda_seq Lambda sequence for regularization
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return List of trained models
#' @export
train_linear_models <- function(X, y, alpha_values = c(0, 0.5, 1),
                                 lambda_seq = NULL, n_folds = 5, seed = 42) {
  
  set.seed(seed)
  
  # Prepare data
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  # Remove zero-variance features
  var_filter <- apply(X, 2, var) > 0
  X <- X[, var_filter, drop = FALSE]
  
  if (ncol(X) == 0) {
    warning("No valid features after variance filtering")
    return(NULL)
  }
  
  models <- list()
  
  # OLS (Ordinary Least Squares)
  tryCatch({
    ols_fit <- lm(y ~ ., data = as.data.frame(X))
    models$ols <- list(
      model = ols_fit,
      type = "ols"
    )
  }, error = function(e) {
    warning("OLS fitting failed: ", e$message)
  })
  
  # Regularized models using glmnet
  for (alpha in alpha_values) {
    model_name <- switch(
      as.character(alpha),
      "0" = "ridge",
      "1" = "lasso",
      "elastic_net"
    )
    
    tryCatch({
      # Fit with cross-validation
      cv_fit <- glmnet::cv.glmnet(
        x = X,
        y = y,
        alpha = alpha,
        nfolds = n_folds,
        lambda = lambda_seq
      )
      
      models[[model_name]] <- list(
        model = cv_fit,
        lambda_min = cv_fit$lambda.min,
        lambda_1se = cv_fit$lambda.1se,
        alpha = alpha,
        type = "glmnet"
      )
    }, error = function(e) {
      warning(sprintf("%s fitting failed: %s", model_name, e$message))
    })
  }
  
  return(models)
}


#' Train Random Forest model
#'
#' @param X Feature matrix (cells × features)
#' @param y Target expression vector
#' @param n_trees Number of trees
#' @param mtry Number of variables to try at each split
#' @param min_node_size Minimum node size
#' @param importance Calculate variable importance
#' @param seed Random seed
#' @return Trained ranger model
#' @export
train_random_forest <- function(X, y, n_trees = 500, mtry = NULL,
                                 min_node_size = 5, importance = "impurity",
                                 seed = 42) {
  
  set.seed(seed)
  
  # Prepare data
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  # Remove zero-variance features
  var_filter <- apply(X, 2, var) > 0
  X <- X[, var_filter, drop = FALSE]
  
  if (ncol(X) == 0) {
    warning("No valid features after variance filtering")
    return(NULL)
  }
  
  # Set mtry if not provided
  if (is.null(mtry)) {
    mtry <- max(1, floor(ncol(X) / 3))
  }
  mtry <- min(mtry, ncol(X))
  
  # Train model
  tryCatch({
    rf_fit <- ranger::ranger(
      x = X,
      y = y,
      num.trees = n_trees,
      mtry = mtry,
      min.node.size = min_node_size,
      importance = importance,
      seed = seed
    )
    
    return(list(
      model = rf_fit,
      type = "ranger",
      mtry = mtry,
      n_trees = n_trees
    ))
  }, error = function(e) {
    warning("Random Forest fitting failed: ", e$message)
    return(NULL)
  })
}


#' Train Deep Neural Network using Keras
#'
#' @param X_train Training features
#' @param y_train Training targets
#' @param X_val Validation features (optional)
#' @param y_val Validation targets (optional)
#' @param hidden_units Vector of hidden layer sizes
#' @param dropout_rate Dropout rate
#' @param learning_rate Learning rate
#' @param batch_size Batch size
#' @param epochs Maximum epochs
#' @param patience Early stopping patience
#' @param seed Random seed
#' @return List with model and training history
#' @export
train_neural_network <- function(X_train, y_train, X_val = NULL, y_val = NULL,
                                  hidden_units = c(256, 128, 64),
                                  dropout_rate = 0.3, learning_rate = 0.001,
                                  batch_size = 256, epochs = 100,
                                  patience = 10, seed = 42) {
  
  # Import keras
  if (!requireNamespace("keras3", quietly = TRUE)) {
    stop("keras3 package required. Install with: install.packages('keras3')")
  }
  
  keras3::set_random_seed(seed)
  
  # Prepare data
  X_train <- as.matrix(X_train)
  y_train <- as.numeric(y_train)
  
  # Remove zero-variance features
  var_filter <- apply(X_train, 2, var) > 0
  X_train <- X_train[, var_filter, drop = FALSE]
  
  if (!is.null(X_val)) {
    X_val <- as.matrix(X_val)[, var_filter, drop = FALSE]
    y_val <- as.numeric(y_val)
  }
  
  n_features <- ncol(X_train)
  
  if (n_features == 0) {
    warning("No valid features after variance filtering")
    return(NULL)
  }
  
  # Build model
  model <- keras3::keras_model_sequential(input_shape = n_features)
  
  for (i in seq_along(hidden_units)) {
    model <- model |>
      keras3::layer_dense(units = hidden_units[i], activation = "relu") |>
      keras3::layer_batch_normalization() |>
      keras3::layer_dropout(rate = dropout_rate)
  }
  
  model <- model |>
    keras3::layer_dense(units = 1, activation = "linear")
  
  # Compile
  model |> keras3::compile(
    optimizer = keras3::optimizer_adam(learning_rate = learning_rate),
    loss = "mse",
    metrics = list("mae")
  )
  
  # Callbacks
  callbacks <- list(
    keras3::callback_early_stopping(
      monitor = if (!is.null(X_val)) "val_loss" else "loss",
      patience = patience,
      restore_best_weights = TRUE
    )
  )
  
  # Train
  validation_data <- if (!is.null(X_val)) list(X_val, y_val) else NULL
  
  history <- model |> keras3::fit(
    x = X_train,
    y = y_train,
    validation_data = validation_data,
    epochs = epochs,
    batch_size = batch_size,
    callbacks = callbacks,
    verbose = 0
  )
  
  return(list(
    model = model,
    history = history,
    type = "keras",
    feature_mask = var_filter
  ))
}


#' Predict using trained model
#'
#' @param model_obj Model object from training functions
#' @param X_new New feature matrix
#' @param lambda Lambda for glmnet models ("min" or "1se")
#' @return Vector of predictions
#' @export
predict_model <- function(model_obj, X_new, lambda = "min") {
  
  X_new <- as.matrix(X_new)
  
  if (model_obj$type == "ols") {
    preds <- predict(model_obj$model, newdata = as.data.frame(X_new))
    
  } else if (model_obj$type == "glmnet") {
    lambda_val <- if (lambda == "min") model_obj$lambda_min else model_obj$lambda_1se
    preds <- as.vector(predict(model_obj$model, newx = X_new, s = lambda_val))
    
  } else if (model_obj$type == "ranger") {
    preds <- predict(model_obj$model, data = X_new)$predictions
    
  } else if (model_obj$type == "keras") {
    # Apply feature mask if needed
    if (!is.null(model_obj$feature_mask)) {
      X_new <- X_new[, model_obj$feature_mask, drop = FALSE]
    }
    preds <- as.vector(model_obj$model |> keras3::predict(X_new, verbose = 0))
    
  } else {
    stop("Unknown model type: ", model_obj$type)
  }
  
  return(preds)
}


#' Calculate prediction metrics
#'
#' @param y_true True values
#' @param y_pred Predicted values
#' @return Named list of metrics
#' @export
calculate_metrics <- function(y_true, y_pred) {
  
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  
  # Remove NA
  valid <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid]
  y_pred <- y_pred[valid]
  
  n <- length(y_true)
  
  # MSE and RMSE
  mse <- mean((y_true - y_pred)^2)
  rmse <- sqrt(mse)
  
  # MAE
  mae <- mean(abs(y_true - y_pred))
  
  # R-squared
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  # Pearson correlation
  if (sd(y_true) > 0 && sd(y_pred) > 0) {
    pearson_r <- cor(y_true, y_pred, method = "pearson")
  } else {
    pearson_r <- NA
  }
  
  # Spearman correlation
  if (length(unique(y_true)) > 1 && length(unique(y_pred)) > 1) {
    spearman_rho <- cor(y_true, y_pred, method = "spearman")
  } else {
    spearman_rho <- NA
  }
  
  return(list(
    n = n,
    mse = mse,
    rmse = rmse,
    mae = mae,
    r2 = r2,
    pearson_r = pearson_r,
    spearman_rho = spearman_rho
  ))
}


#' Train all models for a target gene
#'
#' @param X_train Training features
#' @param y_train Training targets
#' @param X_val Validation features
#' @param y_val Validation targets
#' @param models_to_train Which models to train
#' @param seed Random seed
#' @return List of trained models
#' @export
train_all_models <- function(X_train, y_train, X_val = NULL, y_val = NULL,
                              models_to_train = c("ols", "ridge", "lasso", 
                                                   "elastic_net", "random_forest",
                                                   "neural_network"),
                              seed = 42) {
  
  all_models <- list()
  
  # Linear models
  if (any(c("ols", "ridge", "lasso", "elastic_net") %in% models_to_train)) {
    alpha_vals <- c()
    if ("ridge" %in% models_to_train) alpha_vals <- c(alpha_vals, 0)
    if ("lasso" %in% models_to_train) alpha_vals <- c(alpha_vals, 1)
    if ("elastic_net" %in% models_to_train) alpha_vals <- c(alpha_vals, 0.5)
    
    linear_models <- train_linear_models(
      X = X_train, y = y_train,
      alpha_values = alpha_vals,
      seed = seed
    )
    
    all_models <- c(all_models, linear_models)
  }
  
  # Random Forest
  if ("random_forest" %in% models_to_train) {
    rf_model <- train_random_forest(X = X_train, y = y_train, seed = seed)
    if (!is.null(rf_model)) {
      all_models$random_forest <- rf_model
    }
  }
  
  # Neural Network
  if ("neural_network" %in% models_to_train) {
    nn_model <- train_neural_network(
      X_train = X_train, y_train = y_train,
      X_val = X_val, y_val = y_val,
      seed = seed
    )
    if (!is.null(nn_model)) {
      all_models$neural_network <- nn_model
    }
  }
  
  return(all_models)
}


#' Evaluate all models on test data
#'
#' @param models List of trained models
#' @param X_test Test features
#' @param y_test Test targets
#' @return Data frame of metrics for each model
#' @export
evaluate_all_models <- function(models, X_test, y_test) {
  
  results <- list()
  
  for (model_name in names(models)) {
    model_obj <- models[[model_name]]
    
    tryCatch({
      preds <- predict_model(model_obj, X_test)
      metrics <- calculate_metrics(y_test, preds)
      
      results[[model_name]] <- data.frame(
        model = model_name,
        n_test = metrics$n,
        mse = metrics$mse,
        rmse = metrics$rmse,
        mae = metrics$mae,
        r2 = metrics$r2,
        pearson_r = metrics$pearson_r,
        spearman_rho = metrics$spearman_rho,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Evaluation failed for %s: %s", model_name, e$message))
    })
  }
  
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame())
  }
}


#' Save trained models
#'
#' @param models List of trained models
#' @param output_dir Output directory
#' @param gene_name Target gene name
#' @param modality Modality name
#' @export
save_models <- function(models, output_dir, gene_name, modality) {
  
  model_dir <- file.path(output_dir, modality, "models")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save non-keras models as RDS
  non_keras <- models[sapply(models, function(x) x$type != "keras")]
  if (length(non_keras) > 0) {
    rds_file <- file.path(model_dir, paste0(gene_name, "_models.rds"))
    saveRDS(non_keras, rds_file)
  }
  
  # Save keras models separately
  keras_models <- models[sapply(models, function(x) x$type == "keras")]
  for (name in names(keras_models)) {
    keras_dir <- file.path(model_dir, paste0(gene_name, "_", name, ".keras"))
    keras3::save_model(keras_models[[name]]$model, keras_dir)
  }
}


#' Load trained models
#'
#' @param model_dir Directory containing saved models
#' @param gene_name Target gene name
#' @return List of loaded models
#' @export
load_models <- function(model_dir, gene_name) {
  
  models <- list()
  
  # Load RDS models
  rds_file <- file.path(model_dir, paste0(gene_name, "_models.rds"))
  if (file.exists(rds_file)) {
    models <- c(models, readRDS(rds_file))
  }
  
  # Load keras models
  keras_files <- list.files(model_dir, pattern = paste0(gene_name, "_.*\\.keras$"))
  for (keras_file in keras_files) {
    model_name <- gsub(paste0(gene_name, "_|\\.keras"), "", keras_file)
    keras_path <- file.path(model_dir, keras_file)
    models[[model_name]] <- list(
      model = keras3::load_model(keras_path),
      type = "keras"
    )
  }
  
  return(models)
}
