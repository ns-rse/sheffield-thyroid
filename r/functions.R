#' Function for fitting different models to a dataset
#'
#' For this work we want to fit models to datasets that have been imputed using the MICE package
#'
#' @df data.frame Dataframe to be analysed
#' @imputed logical Whether the data frame is an imputed data set
#' @prop float Proportion to split the data by for training/testing. Default is 0.75
#' @v int Number of folds for cross-validation. Default is 10.
#' @repeats int Number of repeats for cross-validation. Default is 10.
#' @missing_threshold float Proportion of missing data allowed, used in 'recipes::step_filter_missing()'. Default is 0
#' @outcome Outcome variable
#' @continuous vector Vector of continuous predictor variables.
#' @categorical vector Vector of categorical/nominal predictor variables.
#' @model_type str The type of model to fit, options are 'elastic', 'lasso' 'xgboost', 'forest'
#' @mixture float Mixture to use in Elastic Net ('elastic') models. '0' gives Ridge Regression, '1' LASSO. Default
#' '0.5'.
#' @vi_method str The Variable Importance method used to extract importance. Default is 'shap' other options are
#' 'model', 'firm' and 'permute'. See '?vip::vi' for more information.
#' @vi_sort bool Whether to sort the Variable Importance output. See '?vip::vi' for more information. Default FALSE
#' @vi_scale bool Whether to scale the Variable Importance output. See '?vip::vi' for more information. Default FALSE
#' @vi_rank bool Whether to rank the Variable Importance output. See '?vip::vi' for more information. Default TRUE
#' ... Hyperparameters to pass to the selected model_type
#'
#' @seealso [mice, vip]
#' @export
#' @examples
model_fitting <- function(df,
                          imputed = TRUE,
                          prop = 0.75,
                          v = 10,
                          repeats = 10,
                          missing_threshold = 0.0,
                          outcome = "final_pathology",
                          continuous = c(
                              # "albumin",
                              "tsh_value",
                              "lymphocytes",
                              "monocyte",
                              "size_nodule_mm"
                          ),
                          categorical = c(
                              # "ethnicity",
                              "incidental_nodule",
                              "palpable_nodule",
                              "rapid_enlargement",
                              "compressive_symptoms",
                              "hypertension",
                              "vocal_cord_paresis",
                              "graves_disease",
                              "hashimotos_thyroiditis",
                              "family_history_thyroid_cancer",
                              "exposure_radiation",
                              "bta_u_classification",
                              "solitary_nodule",
                              "cervical_lymphadenopathy",
                              "thy_classification"
                          ),
                          model_type = "lasso",
                          mixture = 0.5,
                          vi_method = "model",  ## Ideally want to use 'snap' to make comparisons
                          vi_sort = FALSE,
                          vi_scale = FALSE,
                          vi_rank = TRUE,
                          ...) {
    ## If this is an imputed dataset we need to remove -.imp and -.id variables
    if (imputed) {
        df <- df |>
            dplyr::select(-.imp, -.id)
    }
    ## Select the variables of interest
    df <- df |>
        dplyr::select({{ outcome }}, {{ continuous }}, {{ categorical }})
    ## Split the data into training and test data
    split <- df |>
        rsample::initial_split(prop = prop)
    train <- rsample::training(split)
    test <- rsample::testing(split)
    ## Create k-fold validation on training data
    cv_folds <- rsample::vfold_cv(train, v = v, repeats = repeats)
    ## Setup a recipe, the steps involved removing individuals with missing data
    ## ns-rse 2024-09-20 : Ideally we shouldn't have final_pathology hard coded here but the "embrace" method of
    ## indirection (see vignette("programmin")) doesn't seem to work with recipes. Online example of creating own recipe
    ## step function (https://www.tidymodels.org/learn/develop/recipes/index.html#create-the-function) uses the older
    ## enquos() method
    thyroid_recipe <- recipes::recipe(final_pathology ~ ., data = train) |>
        recipes::step_filter_missing(recipes::all_predictors(), threshold = missing_threshold) |>
        recipes::step_normalize(recipes::all_numeric_predictors()) |>
        recipes::step_dummy(recipes::all_nominal_predictors())
    ## Add recipe to a workrlow
    thyroid_workflow <- workflows::workflow() |>
        workflows::add_recipe(thyroid_recipe)
    ## Conditionally fit the model based on the requested method, this involves setting up a parsnip model specification
    ## and then tune it via tune::tune_grid()
    if (model_type == "elastic") {
        message(paste0("Fitting Elastic Net (mixture : ", mixture, ")", sep=""))
        parsnip_tune <- parsnip::logistic_reg(
            penalty = hardhat::tune(),
            mixture = mixture
        ) |>
            parsnip::set_mode("classification") |>
            parsnip::set_engine("glmnet")
        dials_grid <- dials::grid_regular(dials::penalty(), levels = 50)
    } else if (model_type == "forest") {
        message("Fitting Random Forest")
        ## @ns-rse 2024-09-20 : For flexibility we could set things up to specify different engines here
        parsnip_tune <- parsnip::rand_forest(
            mtry = hardhat::tune(),
            trees = 100,
            min_n = hardhat::tune()
        ) |>
            parsnip::set_mode("classification") |>
            parsnip::set_engine("ranger", importance = "impurity")
        dials_grid <- dials::grid_regular(
            dials::mtry(range = c(5, 10)), # smaller ranges will run quicker
            dials::min_n(range = c(2, 25)),
            levels = 3
        )
    } else if (model_type == "xgboost") {
        message("Fitting Gradient Boosting")
        parsnip_tune <- parsnip::boost_tree(
            mode = "classification",
            trees = 100,
            min_n = hardhat::tune(),
            tree_depth = hardhat::tune(),
            learn_rate = hardhat::tune(),
            loss_reduction = hardhat::tune()
        ) |>
            parsnip::set_mode("classification") |>
            parsnip::set_engine("xgboost", objective = "binary:logistic")
        ## Specify the models tuning parameters using the `dials` package along (https://dials.tidymodels.org/) with the grid
        ## space. This helps identify the hyperparameters with the lowest prediction error.
        dials_params <- dials::parameters(
            dials::min_n(),
            dials::tree_depth(),
            dials::learn_rate(),
            dials::loss_reduction()
        )
        dials_grid <- dials::grid_max_entropy(
            dials_params,
            size = 10
        )
    }
    ## Set common metrics and control
    yardstick_metrics <- yardstick::metric_set(yardstick::roc_auc,
                                               yardstick::accuracy,
                                               yardstick::ppv,
                                               yardstick::npv)
    control <- tune::control_grid(verbose = FALSE)
    ## Tune the model
    tuned_model <- tune::tune_grid(
                             workflows::add_model(thyroid_workflow, parsnip_tune),
                             resamples = cv_folds, ## cv_loo,
                             grid = dials_grid,
                             metrics = yardstick_metrics,
                             control = control
                         )
    ## Select the best model based on ROC Area Under the Curve and finalise the workflow by adding the tuned model and
    ## best metric fit
    best_metric_fit <- tuned_model |>
        tune::select_best()
    final_model <- tune::finalize_workflow(workflows::add_model(thyroid_workflow, spec = parsnip_tune), best_metric_fit)
    fit <- final_model |>
        fit(train)
    ## Extract importance
    importance <- fit |>
        hardhat::extract_fit_parsnip() |>
        vip::vi(method = vi_method,
                sort = vi_sort,
                scale = vi_scale,
                rank = vi_rank)
    ## Add Sign if not present
    if (!("Sign" %in% colnames(importance))) {
        importance <- importance |>
            dplyr::mutate(Sign = "POS")
    }
    ## Plot importance
    importance_plot <- importance |>
        ggplot2::ggplot(mapping = aes(x = Importance, y = Variable, fill = Sign)) +
        ggplot2::geom_col() +
        ggdark::dark_theme_minimal() +
        ggplot2::ggtitle("Train")
    ########################################################
    ## Summarise model on train data set                  ##
    ########################################################
    ## Make classification predictions on train subset
    train[".pred_class"] <- fit |>
        predict(train)
    ## Make probability predictions on train subset
    train[".pred_prob"] <- fit |>
        predict(train, type="prob")
    ## Generate summary statistics from confusion matrix on train subset
    train_summary_metrics <- train |>
        yardstick::conf_mat(final_pathology, .pred_class) |>
        summary()
    ## Extract Area Under the Curve from ROC on train subset
    auc <- train |>
        yardstick::roc_auc(final_pathology, .pred_prob)
    ## Bind summary_metrics and auc into a single data frame
    train_summary_metrics <- rbind(train_summary_metrics, auc)
    ## Extract ROC curve on train subset
    train_roc_curve <- train |>
        ## ToDo - Remove hard coding of final_pathology
        yardstick::roc_curve(final_pathology, .pred_prob)
    ## Plot ROC curve on train subset
    train_roc_curve_plot <- train_roc_curve |>
        ggplot2::ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        ggplot2::geom_path() +
        ggplot2::geom_abline(lty = 3) +
        ggplot2::coord_equal() +
        ggdark::dark_theme_minimal() +
        ggplot2::ggtitle("Train")
    ########################################################
    ## Summarise model on test data set                   ##
    ########################################################
    ## Make classification predictions on test subset
    test[".pred_class"] <- fit |>
        predict(test)
    ## Make probability predictions on test subset
    test[".pred_prob"] <- fit |>
        predict(test, type="prob")
    ## Generate summary statistics from confusion matrix on test subset
    test_summary_metrics <- test |>
        yardstick::conf_mat(final_pathology, .pred_class) |>
        summary()
    ## Extract Area Under the Curve from ROC on test subset
    auc <- test |>
        yardstick::roc_auc(final_pathology, .pred_prob)
    ## Bind summary_metrics and auc into a single data frame
    test_summary_metrics <- rbind(test_summary_metrics, auc)
    ## Extract ROC curve on test subset
    test_roc_curve <- test |>
        ## ToDo - Remove hard coding of final_pathology
        yardstick::roc_curve(final_pathology, .pred_prob)
    ## Plot ROC curve on test subset
    test_roc_curve_plot <- test_roc_curve |>
        ggplot2::ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        ggplot2::geom_path() +
        ggplot2::geom_abline(lty = 3) +
        ggplot2::coord_equal() +
        ggdark::dark_theme_minimal() +
        ggplot2::ggtitle("Test")
    ## Combine results into a single list for returning
    results <- list()
    results$train <- train
    results$test <- test
    results$folds <- cv_folds
    results$recipe <- recipe
    results$workflow <- workflow
    results$parsnip <- parsnip_tune
    results$grid <- dials_grid
    results$yardstick <- yardstick_metrics
    results$control <- control
    results$best <- best_metric_fit
    results$final_model <- final_model
    results$importance <- importance
    results$importance_plot <- importance_plot
    if (model_type == "elastic") {
        results$broom_estimate <- fit |> broom::tidy()
    }
    results$fit <- fit
    results$train_summary_metrics <- train_summary_metrics
    results$train_roc_curve <- train_roc_curve
    results$train_roc_curve_plot <- train_roc_curve_plot
    results$test_summary_metrics <- test_summary_metrics
    results$test_roc_curve <- test_roc_curve
    results$test_roc_curve_plot <- test_roc_curve_plot
    results
}


#' Extract and combine summary statistics from analysis of multiple imputed data sets
#'
#' @imputed_models list A list of imputed models from having mapped model_fitting() with purrr() on an imputed data set.
#' @model_type str The type of model that has been fitted across the imputed datasets
#'
tidy_imputed_models <- function(imputed_models, model_type) {
    results <- list()
    ## Extract importance
    importance <- imputed_models |>
        purrr::map(`[[`, "importance") |>
        dplyr::bind_rows(.id = "imputation") |>
        ## tidyr::spread(key = imputation, value = Importance)
        tidyr::pivot_wider(names_from = imputation, values_from = Importance)
    if (model_type %in% c("elastic", "xgboost")){
        importance <- importance |>
            dplyr::arrange(Variable, Sign)
    }
    else {
        importance <- importance |>
            dplyr::arrange(Variable)
    }
    ## Extract Summary Statistics for test and train
    test_metrics <- imputed_models |>
        purrr::map(`[[`, "test_summary_metrics") |>
        dplyr::bind_rows(.id = "imputation") |>
        ## dplyr::select(-.estimate) |>
        ## tidyr::spread(key = imputation, value = .estimator)
        tidyr::pivot_wider(names_from = imputation, values_from = .estimate)
    train_metrics <- imputed_models |>
        purrr::map(`[[`, "train_summary_metrics") |>
        dplyr::bind_rows(.id = "imputation") |>
        ## dplyr::select(-.estimate) |>
        tidyr::spread(key = imputation, value = .estimate)
    ## Extract ROC curves (NB these are not converted to wide format to simplify plotting)
    test_roc_metrics <- imputed_models |>
        purrr::map(`[[`, "test_roc_curve") |>
        dplyr::bind_rows(.id = "imputation")
    train_roc_metrics <- imputed_models |>
        purrr::map(`[[`, "train_roc_curve") |>
        dplyr::bind_rows(.id = "imputation")
    ## TODO - Plot test and train ROC curves on the same graph coloured by imputation
    ## Add all summaries to results
    results$importance <- importance
    results$test_metrics <- test_metrics
    results$train_metrics <- train_metrics
    results$test_roc_metrics <- test_roc_metrics
    results$train_roc_metrics <- train_roc_metrics
    results
}
