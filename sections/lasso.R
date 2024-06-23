## ----------------------------------------------------------------------------------------
#| label: lasso
#| purl: true
#| eval: true
#| echo: true
#| output: false
#| cache: true
## Specify the LASSO model using parsnip, the key here is the use of the glmnet engine which is the R package for
## fitting LASSO regression. Technically the package fits Elastic Net but with a mixture value of 1 it is equivalent to
## a plain LASSO (mixture value of 0 is equivalent to Ridge Regression in an Elastic Net)
tune_spec_lasso <- parsnip::logistic_reg(penalty = hardhat::tune(), mixture = 1) |>
  parsnip::set_engine("glmnet")

## Tune the LASSO parameters via cross-validation
lasso_grid <- tune::tune_grid(
  object = workflows::add_model(thyroid_workflow, tune_spec_lasso),
  resamples = cv_folds,
  grid = dials::grid_regular(penalty(), levels = 50)
)


## ----------------------------------------------------------------------------------------
#| label: fig-lasso-tune-autoplot
#| fig-cap: Autoplot of LASSO grid search
#| purl: true
#| eval: true
#| echo: true
#| output: true
## Plot the tuning search results, see https://tune.tidymodels.org/reference/autoplot.tune_results.html
##
## This shows how the evaluation metrics change over time with each iteration of the Lasso (we know we are running a
## LASSO because when we setup tune_spec_lasso the mixture = 1 when we define parsnip::logistic_reg()
tune::autoplot(lasso_grid)


## ----------------------------------------------------------------------------------------
#| label: fig-lasso-kfold-eval-importance
#| fig-cap: Importance of variables fitted using LASSO
#| purl: true
#| eval: true
#| echo: true
#| output: true
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


## ----------------------------------------------------------------------------------------
#| label: lasso-save
#| purl: true
#| eval: true
#| echo: true
#| output: false
save(tune_spec_lasso, lasso_grid, final_lasso_kfold, lasso_kfold_roc_auc,
  file = "data/r/lasso.RData"
)

