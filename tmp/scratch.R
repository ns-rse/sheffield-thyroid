library(dplyr)
library(fastshap)
library(ggplot2)
library(ggdark)
library(gtsummary)
library(Hmisc)
library(knitr)
library(readr)
library(rmarkdown)
library(visdat)
library(naniar)
library(purrr)

## Libraries for Tidymodelling
library(dials)
library(kernlab)
library(knitr)
library(tidymodels)
library(tidyverse)
library(vip)

base_dir <- getwd()
data_dir <- paste(base_dir, "data", sep="/")
r_dir <- paste(data_dir, "r", sep="/")
## Load and retain only those with final_pathology
df <- readRDS(paste(r_dir, "clean.rds", sep = "/"))
df_complete <- readRDS(paste(r_dir, "df_complete.rds", sep = "/"))
mice_pmm <- readRDS(file = paste(r_dir, "mice_pmm.rds", sep = "/"))
mice_cart <- readRDS(file = paste(r_dir, "mice_cart.rds", sep = "/"))
mice_rf <- readRDS(file = paste(r_dir, "mice_rf.rds", sep = "/"))
thyroid_workflow <- readRDS(file = paste(r_dir, "thyroid_workflow.rds", sep = "/"))
####################################################
## 2024-09-20 - Purrr alternative, write a function
##              that does the modelling and map that
####################################################
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
#' @model_type str The type of model to fit, options are 'elastic', 'lasso' 'xgboost', 'svm'
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
                          continuous = c("albumin", "tsh_value", "lymphocytes", "monocyte", "size_nodule_mm"),
                          categorical = c(
                              "ethnicity",
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
    ## @ns-rse 2024-09-20 : Refactor and abstract out common components have parsnip::() and grid_regular() set in each
    ## conditional section and a final tune::tune_grid() that takes whatever these values are set to.
    if (model_type == "elastic") {
        message(paste0("Fitting Elastic Net (mixture : ", mixture, ")", sep=""))
        parsnip_tune <- parsnip::logistic_reg(
            penalty = hardhat::tune(),
            mixture = mixture
        ) |>
            parsnip::set_mode("classification") |>
            parsnip::set_engine("glmnet")
        dials_grid <- dials::grid_regular(dials::penalty(), levels = 50)
        ## grid <- tune::tune_grid(
        ##     object = workflows::add_model(thyroid_workflow, parsnip_tune),
        ##     resamples = cv_folds,
        ##     grid = dials_grid
        ## )
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
        ## grid <- tune::tune_grid(
        ##     object = workflows::add_model(thyroid_workflow, parsnip_tune),
        ##     resamples = cv_folds,
        ##     grid = dials_grid
        ## )
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
    ## Not bothering with SVM modelling, there are errors in the calculations during modelling and no methods available
    ## to extract Variable Importance
    ## else if (model_type == "svm") {
    ##     message("Fitting Support Vector Machine")
    ##     parsnip_tune <- parsnip::svm_rbf(cost = tune()) |>
    ##       parsnip::set_mode("classification") |>
    ##       parsnip::set_engine("kernlab")
    ##     dials_grid <- dials::grid_regular(dials::cost(), levels = 20)
    ## }
    ## Set common metrics and control
    yardstick_metrics <- yardstick::metric_set(yardstick::roc_auc,
                                               yardstick::accuracy,
                                               yardstick::ppv)
    control <- tune::control_grid(verbose = FALSE)
    ## Tune the model
    tuned_model <- tune::tune_grid(
                             workflows::add_model(thyroid_workflow, parsnip_tune),
                             resamples = cv_folds, ## cv_loo,
                             grid = dials_grid,
                             metrics = yardstick_metrics,
                             control = control
                         )
    ## Select the best model based on ROC Area Under the Curve and finalise the workflow by adding the tunde model and
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
        ggdark::dark_theme_minimal()
    ########################################################
    ## Summarise model on train data set                  ##
    ########################################################
    ## Make predictions on train subset
    train[".pred_class"] <- fit |>
        predict(train)
    train[".pred_prob"] <- fit |>
        predict(train, type="prob")
    ## ToDo - Remove hard coding of final_pathology from yardstick calls (enquotes?)
    train_summary_metrics <- train |>
        yardstick::conf_mat(final_pathology, .pred_class) |>
        summary()
    auc <- train |>
        yardstick::roc_auc(final_pathology, .pred_prob)
    train_summary_metrics <- rbind(train_summary_metrics, auc)
    train_roc_curve <- train |>
        ## ToDo - Remove hard coding of final_pathology
        yardstick::roc_curve(final_pathology, .pred_prob)
    train_roc_curve_plot <- train_roc_curve |>
        ggplot2::ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        ggplot2::geom_path() +
        ggplot2::geom_abline(lty = 3) +
        ggplot2::coord_equal() +
        ggdark::dark_theme_minimal()
    ########################################################
    ## Summarise model on train data set                  ##
    ########################################################
    ## Make predictions on test subset
    test[".pred_class"] <- fit |>
        predict(test)
    test[".pred_prob"] <- fit |>
        predict(test, type="prob")
    ## ToDo - Remove hard coding of final_pathology from yardstick calls (enquotes?)
    test_summary_metrics <- test |>
        yardstick::conf_mat(final_pathology, .pred_class) |>
        summary()
    auc <- test |>
        yardstick::roc_auc(final_pathology, .pred_prob)
    test_summary_metrics <- rbind(test_summary_metrics, auc)
    test_roc_curve <- test |>
        ## ToDo - Remove hard coding of final_pathology
        yardstick::roc_curve(final_pathology, .pred_prob)
    test_roc_curve_plot <- test_roc_curve |>
        ggplot2::ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        ggplot2::geom_path() +
        ggplot2::geom_abline(lty = 3) +
        ggplot2::coord_equal() +
        ggdark::dark_theme_minimal()
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
elastic_model <- mice_cart$imputed |>
  dplyr::filter(.imp == 1) |>
  ## model_fitting(model_type = "elastic", mixture = 1)   ## LASSO
  model_fitting(model_type = "elastic", mixture = 0.5) ## Elastic Net

rforest_model <- mice_cart$imputed |>
  dplyr::filter(.imp == 1) |>
  model_fitting(model_type = "forest")

lasso_model <- mice_cart$imputed |>
  dplyr::filter(.imp == 1) |>
  model_fitting(model_type = "elastic", mixture = 1)   ## LASSO

xgb_model <- mice_cart$imputed |>
  dplyr::filter(.imp == 1) |>
  model_fitting(model_type = "xgboost")


## Map across imputed data sets.
## 1. Extract the data from the mice object
## 2. Remove Original data set.
## 3. Remove the "Original" value from .imp levels
## 4. Save to cart_imputed
cart_imputed <- mice_cart$imputed |>
    dplyr::filter(.imp != "Original") |>
    dplyr::mutate(.imp = droplevels(.imp))

## 1. Split the cart_imputed data by the .imp variable (which indicates which imputation the data is from)
## 2. The first argument to purrr::map(), is an imputed data set and referenced by '\(df)' this means we can use 'df' in
##    the function we are mapping to, in this case model_fitting() so it will receive an imputed dataset.
## 3. The rest of the arguments to model_fitting() follow.
imputed_cart_model_elastic <- cart_imputed |>
    split(cart_imputed$.imp) |>
    purrr::map(\(df) model_fitting(df,
        model_type = "elastic",
        mixture = 0.5,
        imputed = TRUE
    ))


## Combine and tabulate (very hacky, don't have time to write elegant code)...
## - importance
## - summary statistics
## - ROC
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
    test_metrics <- imputed_cart_model |>
        purrr::map(`[[`, "test_summary_metrics") |>
        dplyr::bind_rows(.id = "imputation") |>
        ## dplyr::select(-.estimate) |>
        ## tidyr::spread(key = imputation, value = .estimator)
        tidyr::pivot_wider(names_from = imputation, values_from = .estimate)
    train_metrics <- imputed_cart_model |>
        purrr::map(`[[`, "train_summary_metrics") |>
        dplyr::bind_rows(.id = "imputation") |>
        ## dplyr::select(-.estimate) |>
        tidyr::spread(key = imputation, value = .estimate)
    ## Extract ROC curves (NB these are not converted to wide format to simplify plotting)
    test_roc_metrics <- imputed_cart_model |>
        purrr::map(`[[`, "test_roc_curve") |>
        dplyr::bind_rows(.id = "imputation")
    train_roc_metrics <- imputed_cart_model |>
        purrr::map(`[[`, "train_roc_curve") |>
        dplyr::bind_rows(.id = "imputation")
    ## TODO - Plot test and train ROC curves coloured by imputation

    ## Add all summaries to results
    results$importance <- importance
    results$test_metrics <- test_metrics
    results$train_metrics <- train_metrics
    results$test_roc_metrics <- test_roc_metrics
    results$train_roc_metrics <- train_roc_metrics
    results
}
test <- tidy_imputed_models(imputed_cart_model_elastic)
test$test_roc_metrics |>
    ## dplyr::filter(imputation == 1) |>
    dplyr::arrange(imputation, .threshold, -specificity) |>
    ggplot2::ggplot(aes(1 - specificity, sensitivity, color = imputation)) + geom_line()

## ToDo...
##
##
## XGBoost : Need to get feature names from importance if we want to use Shap values/method (useful for making
## comparisons of importance).
## Combine estimates by averaging using the AIC as weight, need to extract likelihood for each model may be able to
## leverage mitools (https://cran.r-project.org/web/packages/mitools/index.html) although we could potentially use
## mice::pool() to get the estimates directly and then mitools to get p-values see
## https://francish.net/post/multiple-imputation-in-r-with-regression-output/

####################################################
## 2024-09-19 - Purrr with MICE
####################################################
library(tidyverse)
library(tidymodels)
library(mice)
## Make some missing data
iris[16, 2] <- as.numeric(NA)
iris[46, 2] <- as.numeric(NA)
iris[76, 2] <- as.numeric(NA)

iris[96, 2] <- as.numeric(NA)
iris[106, 2] <- as.numeric(NA)
iris[18, 3] <- as.numeric(NA)
iris[48, 3] <- as.numeric(NA)
iris[78, 3] <- as.numeric(NA)
iris[98, 3] <- as.numeric(NA)
iris[108, 3] <- as.numeric(NA)
iris[9, 4] <- as.numeric(NA)
iris[29, 4] <- as.numeric(NA)
iris[69, 4] <- as.numeric(NA)
iris[89, 4] <- as.numeric(NA)
iris[109, 4] <- as.numeric(NA)
iris[129, 4] <- as.numeric(NA)

## Impute missing values
imputations <- 4
iterations <- 5
seed <- 238992
mids <- iris |>
  dplyr::select(-Sepal.Length) |>
  mice::futuremice(
    m = imputations,
    method = "cart",
    iterations = iterations,
    n.core = imputations,
    parallelseed = seed
  )
## Generate imputed dataset
imputed <- mids |>
  mice::complete(action = "long", include = TRUE)
## Convert .imp variable which indicates the imputation set to factor with original dataset labelled as such
imputed <- imputed |>
  dplyr::mutate(.imp = factor(.imp,
    levels = seq(0, imputations, 1),
    labels = c("Original", as.character(seq(1, imputations, 1)))
  ))
## Bind the Sepal.length to each dataset of imputed (including original)
outcome = iris[["Sepal.Length"]]
n <- imputations
while (n > 0) {
  outcome <- append(outcome, iris[["Sepal.Length"]])
  n <- n - 1
}
imputed <- cbind(imputed, outcome)
## Sort out column names
colnames(imputed) <- stringr::str_replace(colnames(imputed), "outcome", "Sepal.Length")

## We now have six data sets in imputed differentiated by the .imp variable.names

## Make a recipe and workflow to analyse each of these.
set.seed(5039378)
## Split the data
purrr_split_iris <- imputed |>
    dplyr::select(-.id) |>
    split(imputed$.imp) |>
    purrr::map(\(df) rsample::initial_split(df, prop = 0.75))
purrr_train_iris <- purrr_split_iris |>
    purrr::map(\(split) rsample::training(split))
purrr_test_iris <- purrr_split_iris |>
    purrr::map(\(split) rsample::testing(split))

## Define k-fold cross validation
purrr_cv_iris <- purrr_train_iris|>
    purrr::map(\(df) rsample::vfold_cv(df, v = 10, repeats = 10))

## Define a recipe adding steps to normalize numeric variables and create dummies for nominal variables
purrr_recipe_iris <- purrr_train_iris |>
    purrr::map(\(df) recipes::recipe(Sepal.Length ~ ., data = df)) |>
    purrr::map(\(df) recipes::step_normalize(df, recipes::all_numeric_predictors())) |>
    purrr::map(\(df) recipes::step_dummy(df, recipes::all_nominal_predictors()))

####################################################
## 2024-09-03 - Purrr with MICE
####################################################
## Load imputed data set (two imputations)
mice_cart <- readRDS(paste(r_dir, "mice_cart.rds", sep = "/"))
purrr_split <- mice_cart$imputed |>
    split(mice_cart$imputed$.imp) |>
    purrr::map(\(df) rsample::initial_split(df, prop = 0.75))

purrr_train <- purrr_split |>
    purrr::map(\(split) rsample::training(split))
purrr_test <- purrr_split |>
    purrr::map(\(split) rsample::testing(split))

split <- mice_rf$imputed |>
    dplyr::filter(.imp == 1) |>
    dplyr::select(-.imp, -.id) |>
    rsample::initial_split(prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)

####################################################
## 2024-07-19 - Why aren't different imputation methods giving different numbers
####################################################

## Old method
test_pmm <- df_complete |>
    dplyr::select(-final_pathology) |>    ## We deliberately exclude the outcome of interest
    mice::futuremice(m = imputations, method="pmm", n.core=cores) |>
    mice::complete(action = "long", include = TRUE)
test_pmm <- dplyr::mutate(test_pmm, .imp = factor(.imp, levels = c(0,1,2,3,4), labels = c("Original", "1", "2", "3", "4")))
## CART (Classification And Regression Trees)
test_cart <- df_complete |>
    dplyr::select(-final_pathology) |>    ## We deliberately exclude the outcome of interest
    mice::futuremice(m = imputations, method="cart", n.core=cores) |>
    mice::complete(action = "long", include = TRUE)
test_cart <- dplyr::mutate(test_cart, .imp = factor(.imp, levels = c(0,1,2,3,4), labels = c("Original", "1", "2", "3", "4")))
## Random Forest
test_rf <- df_complete |>
    dplyr::select(-final_pathology) |>    ## We deliberately exclude the outcome of interest
    mice::futuremice(m = imputations, method="rf", n.core=cores) |>
    mice::complete(action = "long", include = TRUE)
test_rf <- dplyr::mutate(test_rf, .imp = factor(.imp, levels = c(0,1,2,3,4), labels = c("Original", "1", "2", "3", "4")))

####################################################
## 2024-07-18 - Investigating with recipes::step_unknown()
####################################################
split <- df_complete|>
    rsample::initial_split(prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)

## Without removing anything what doe we have
thyroid_recipe <- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
str(thyroid_recipe)

## Setup a recipe and workflow that...
##    + normalize numeric variables
##    + creates dummy variables for nominal (binary/continuous) variables
##    + removes variables that have correlation > 0.9
thyroid_recipe_step_unknown <- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_unknown(recipes::all_nominal_predictors(), new_level = "unknown") |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
str(thyroid_recipe_step_unknown)

## Compare the packages
waldo::compare(thyroid_recipe, thyroid_recipe_step_unknown)

## Try fitting a model
thyroid_recipe_step_unknown <- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_unknown(recipes::all_nominal_predictors(), new_level = "unknown") |>
    recipes::step_dummy(recipes::all_nominal_predictors())
thyroid_workflow <- workflows::workflow() |>
  workflows::add_recipe(thyroid_recipe_step_unknown)
cv_folds <- rsample::vfold_cv(train, v = 10, repeats = 10)
tune_spec_lasso <- parsnip::logistic_reg(penalty = hardhat::tune(), mixture = 1) |>
  parsnip::set_engine("glmnet")
## Tune the LASSO parameters via cross-validation
lasso_grid <- tune::tune_grid(
  object = workflows::add_model(thyroid_workflow, tune_spec_lasso),
  resamples = cv_folds,
  grid = dials::grid_regular(penalty(), levels = 50)
)

####################################################
## 2024-07-18 - Investigating with recipes::step_corr()
####################################################
split <- df_complete|>
    rsample::initial_split(prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)

## Without removing anything what doe we have
thyroid_recipe <- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
str(thyroid_recipe)

## Setup a recipe and workflow that...
##    + normalize numeric variables
##    + creates dummy variables for nominal (binary/continuous) variables
##    + removes variables that have correlation > 0.9
thyroid_recipe_step_corr <- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_corr(all_predictors(), threshold = 0.9) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)
str(thyroid_recipe_step_corr)

## Compare the packages
waldo::compare(thyroid_recipe, thyroid_recipe_step_corr)

####################################################
## 2024-07-18 - Investigating with recipes::prep() and recipes::bake()
####################################################
split <- df_complete|>
    rsample::initial_split(prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)

## Setup a recipe and workflow with all steps
thyroid_recipe_all<- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep()
all <- thyroid_recipe_all |> bake(new_data = NULL)
str(all)

## Setup a recipe and workflow with all steps increasing threshold for missing predictors
thyroid_recipe_high_missing_threshold<- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_filter_missing(recipes::all_predictors(), threshold = 0.5) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep()
high_missing_threshold <- thyroid_recipe_high_missing_threshold |> bake(new_data = NULL)
str(high_missing_threshold)

## Setup a recipe and workflow with all steps increasing threshold for missing predictors
thyroid_recipe_skip_filter_missing <- recipes::recipe(final_pathology ~ ., data = train) |>
    ## recipes::step_filter_missing(recipes::all_predictors(), threshold = 0, skip = TRUE) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep()
skip_filter_missing<- thyroid_recipe_skip_filter_missing |> bake(new_data = NULL)
str(skip_filter_missing)

## Setup a recipe and workflow with all steps increasing threshold for missing predictors
thyroid_recipe_investigate_missing_threshold<- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_filter_missing(recipes::all_predictors(), threshold = 0.1) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep()
investigate_missing_threshold <- thyroid_recipe_investigate_missing_threshold |> bake(new_data = NULL)
str(investigate_missing_threshold)

## Test fitting a model using higher threshold
thyroid_recipe_investigate_missing_threshold<- recipes::recipe(final_pathology ~ ., data = train) |>
    recipes::step_filter_missing(recipes::all_predictors(), threshold = 0.1) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors())
thyroid_workflow <- workflows::workflow() |>
  workflows::add_recipe(thyroid_recipe_investigate_missing_threshold)
cv_folds <- rsample::vfold_cv(train, v = 10, repeats = 10)
tune_spec_lasso <- parsnip::logistic_reg(penalty = hardhat::tune(), mixture = 1) |>
  parsnip::set_engine("glmnet")
## Tune the LASSO parameters via cross-validation
lasso_grid <- tune::tune_grid(
  object = workflows::add_model(thyroid_workflow, tune_spec_lasso),
  resamples = cv_folds,
  grid = dials::grid_regular(penalty(), levels = 50)
)
## K-fold best fit for LASSO
lasso_kfold_roc_auc <- lasso_grid |>
  tune::select_best(metric = "roc_auc")


## Fit the final LASSO model
final_lasso_kfold <- tune::finalize_workflow(
  workflows::add_model(thyroid_workflow, tune_spec_lasso),
  lasso_kfold_roc_auc)

final_lasso_kfold |>
  fit(train) |>
  hardhat::extract_fit_parsnip() |>
  vip::vi(lambda = lasso_kfold_roc_auc$penalty) |>
  dplyr::mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) |>
  ggplot(mapping = aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  dark_theme_minimal()
ggplot2::ggsave("tmp/lasso_exclude_step_filter_missing.png")

####################################################
## 2024-07-18 - Testing without the `step_filter_missing()`
####################################################
split <- df_complete|>
    rsample::initial_split(prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)

## Setup cross-validation with k-folds (n = 10)
cv_folds <- rsample::vfold_cv(train, v = 10, repeats = 10)


## Setup a recipe and workflow
thyroid_recipe <- recipes::recipe(final_pathology ~ ., data = train) |>
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::prep()
thyroid_workflow <- workflows::workflow() |>
  workflows::add_recipe(thyroid_recipe)

## Try tuning a LASSO model
tune_spec_lasso <- parsnip::logistic_reg(penalty = hardhat::tune(), mixture = 1) |>
  parsnip::set_engine("glmnet")

## Tune the LASSO parameters via cross-validation
lasso_grid <- tune::tune_grid(
  object = workflows::add_model(thyroid_workflow, tune_spec_lasso),
  resamples = cv_folds,
  grid = dials::grid_regular(penalty(), levels = 50)
)

## K-fold best fit for LASSO
lasso_kfold_roc_auc <- lasso_grid |>
  tune::select_best(metric = "roc_auc")


## ----------------------------------------------------------------------------------------------------------------------------------
## Fit the final LASSO model
final_lasso_kfold <- tune::finalize_workflow(
  workflows::add_model(thyroid_workflow, tune_spec_lasso),
  lasso_kfold_roc_auc)

final_lasso_kfold |>
  fit(train) |>
  hardhat::extract_fit_parsnip() |>
  vip::vi(lambda = lasso_kfold_roc_auc$penalty) |>
  dplyr::mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) |>
  ggplot(mapping = aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  dark_theme_minimal()
ggplot2::ggsave("tmp/lasso_exclude_step_filter_missing.png")

####################################################
## 2024-07-18 - Testing with an imputed data set
####################################################
##
## Impute data
imputations = 1
cores = imputations
mice_pmm_single <- df_complete |>
    dplyr::select(-final_pathology) |> ## We deliberately exclude the outcome of interest
    mice::futuremice(m = imputations, method = "pmm", n.core = cores) |>
    mice::complete(action = "long", include = FALSE)
## Bind the final_pathology, removing the .id and .imp columns too and renaming columns
mice_pmm_single <- cbind(mice_pmm_single, df_complete$final_pathology) |>
    dplyr::select(-.id)
colnames(mice_pmm_single) <- stringr::str_replace(colnames(mice_pmm_single), "df_complete\\$", "")

## Split a single imputed data set into train/test
split <- mice_pmm_single |>
    dplyr::filter(.imp == 1) |>
    dplyr::select(-.imp) |>
    rsample::initial_split(prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)

## Setup cross-validation with k-folds (n = 10)
cv_folds <- rsample::vfold_cv(train, v = 10, repeats = 10)


## Setup a recipe and workflow
thyroid_recipe <- recipes::recipe(final_pathology ~ ., data = train) |>
  ## @ns-rse 2024-06-14 :
  ## This step can be used to filter observations with missing data, see the manual pages for more details
  ## https://recipes.tidymodels.org/reference/step_filter_missing.html
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  ## @ns-rse 2024-06-14 :
  ## We first normalise the data _before_ we generate dummies otherwise the dummies, which are numerical, get normalised
  ## too
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
  recipes::step_dummy(recipes::all_nominal_predictors())
thyroid_workflow <- workflows::workflow() |>
  workflows::add_recipe(thyroid_recipe)

## Try tuning a LASSO model
tune_spec_lasso <- parsnip::logistic_reg(penalty = hardhat::tune(), mixture = 1) |>
  parsnip::set_engine("glmnet")

## Tune the LASSO parameters via cross-validation
lasso_grid <- tune::tune_grid(
  object = workflows::add_model(thyroid_workflow, tune_spec_lasso),
  resamples = cv_folds,
  grid = dials::grid_regular(penalty(), levels = 50)
)

## K-fold best fit for LASSO
lasso_kfold_roc_auc <- lasso_grid |>
  tune::select_best(metric = "roc_auc")


## ----------------------------------------------------------------------------------------------------------------------------------
## Fit the final LASSO model
final_lasso_kfold <- tune::finalize_workflow(
  workflows::add_model(thyroid_workflow, tune_spec_lasso),
  lasso_kfold_roc_auc)

final_lasso_kfold |>
  fit(train) |>
  hardhat::extract_fit_parsnip() |>
  vip::vi(lambda = lasso_kfold_roc_auc$penalty) |>
  dplyr::mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) |>
  ggplot(mapping = aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  dark_theme_minimal()
ggplot2::ggsave("tmp/lasso_from_imputed.png")

##########################################

## 2024-07-11 - Investigating why Tidymodels aren't using all data
## Split data into train/test
split <- rsample::initial_split(df_complete, prop = 0.75)
train <- rsample::training(split)
test <- rsample::testing(split)


## Setup cross-validation with k-folds (n = 10)
cv_folds <- rsample::vfold_cv(train, v = 10, repeats = 10)

## Run glm() on subsets of variables
##
## Clinical variables
glm_clin <- glm(final_pathology ~ age_at_scan + gender + incidental_nodule + palpable_nodule +
    rapid_enlargement + compressive_symptoms + hashimotos_thyroiditis + family_history_thyroid_cancer + exposure_radiation +
    tsh_value + size_nodule_mm + solitary_nodule + cervical_lymphadenopathy,
    data = train,
    family = binomial(link = "logit"))
summary(glm_clin)
glm_all

## Biochemical markers
glm_biochem <- glm(final_pathology ~ age_at_scan + gender + tsh_value + albumin + lymphocytes + monocyte,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_biochem)
glm_biochem

## Ultrasound features
glm_ultrasound <- glm(
    final_pathology ~ age_at_scan + gender + incidental_nodule + tsh_value + size_nodule_mm +
        solitary_nodule + cervical_lymphadenopathy + bta_u_classification, # + thy_classification,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_ultrasound)
glm_ultrasound

## Saturated model
glm_all <- glm(
    final_pathology ~ age_at_scan +
        gender +
        incidental_nodule +
        palpable_nodule +
        rapid_enlargement +
        compressive_symptoms +
        hashimotos_thyroiditis +
        family_history_thyroid_cancer +
        ## exposure_radiation +
        ## tsh_value +
        size_nodule_mm +
        solitary_nodule +
        cervical_lymphadenopathy +
        ## albumin +
        ## lymphocytes +
        ## monocyte +
        incidental_nodule +
        size_nodule_mm +
        solitary_nodule +
        bta_u_classification,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_all)
glm_all


## Just incidental and palpable nodules
glm_incidental_palpable <- glm(
    final_pathology ~
        incidental_nodule +
        palpable_nodule,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_incidental_palpable)
glm_incidental_palpable


####################################################
## Lets try out the tidy modelling                ##
####################################################

## incidental nodule
glm_incidental <- glm(
    final_pathology ~ incidental_nodule,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_incidental)
glm_incidental

recipe_incidental <- recipes::recipe(final_pathology ~ incidental_nodule, data = train) |>
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
  recipes::step_dummy(recipes::all_nominal_predictors())
## Workflow
workflow_incidental <- workflows::workflow() |>
  workflows::add_recipe(recipe_incidental)
## Model
logistic_model <- logistic_reg() |>
  set_engine("glm")
log_incidental <- workflow_incidental |>
  add_model(logistic_model)
## Fit
log_incidental_fit <- fit(log_incidental, data = train)
log_incidental_fit


## palpable nodule
glm_palpable <- glm(
    final_pathology ~ palpable_nodule,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_palpable)
glm_palpable


recipe_palpable <- recipes::recipe(final_pathology ~ palpable_nodule, data = train) |>
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
  recipes::step_dummy(recipes::all_nominal_predictors())
## Workflow
workflow_palpable <- workflows::workflow() |>
  workflows::add_recipe(recipe_palpable)
## Model
logistic_model <- logistic_reg() |>
  set_engine("glm")
log_palpable <- workflow_palpable |>
  add_model(logistic_model)
## Fit
log_palpable_fit <- fit(log_palpable, data = train)
log_palpable_fit

## incidental nodule and palpable nodule
glm_incidental_palpable <- glm(
    final_pathology ~ incidental + palpable_nodule,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_incidental_palpable)
glm_incidental_palpable


recipe_incidental_palpable <- recipes::recipe(final_pathology ~ incidental_nodule + palpable_nodule, data = train) |>
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
  recipes::step_dummy(recipes::all_nominal_predictors())
## Workflow
workflow_incidental_palpable <- workflows::workflow() |>
  workflows::add_recipe(recipe_incidental_palpable)
## Model
logistic_model <- logistic_reg() |>
  set_engine("glm")
log_incidental_palpable <- workflow_incidental_palpable |>
  add_model(logistic_model)
## Fit
log_incidental_palpable_fit <- fit(log_incidental_palpable, data = train)
log_incidental_palpable_fit



## incidental nodule and age nodule
glm_age_incidental <- glm(
    final_pathology ~ age_at_scan+ incidental_nodule,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_age_incidental)
glm_age_incidental


recipe_age_incidental <- recipes::recipe(final_pathology ~ age_at_scan + incidental_nodule, data = train) |>
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
  recipes::step_dummy(recipes::all_nominal_predictors())
## Workflow
workflow_age_incidental <- workflows::workflow() |>
  workflows::add_recipe(recipe_age_incidental)
## Model
logistic_model <- logistic_reg() |>
  set_engine("glm")
log_age_incidental <- workflow_age_incidental |>
  add_model(logistic_model)
## Fit
log_age_incidental_fit <- fit(log_age_incidental, data = train)
log_age_incidental_fit



## palpable nodule and age nodule
glm_age_palpable <- glm(
    final_pathology ~ age_at_scan+ palpable_nodule,
    data = train,
    family = binomial(link = "logit")
)
summary(glm_age_palpable)
glm_age_palpable


recipe_age_palpable <- recipes::recipe(final_pathology ~ age_at_scan + palpable_nodule, data = train) |>
  recipes::step_filter_missing(recipes::all_predictors(), threshold = 0) |>
  recipes::step_normalize(recipes::all_numeric_predictors()) |>
  recipes::step_dummy(recipes::all_nominal_predictors())
## Workflow
workflow_age_palpable <- workflows::workflow() |>
  workflows::add_recipe(recipe_age_palpable)
## Model
logistic_model <- logistic_reg() |>
  set_engine("glm")
log_age_palpable <- workflow_age_palpable |>
  add_model(logistic_model)
## Fit
log_age_palpable_fit <- fit(log_age_palpable, data = train)
log_age_palpable_fit
