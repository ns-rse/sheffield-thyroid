# Sheffield Thyroid Cancer Study

This repository contains the [Quarto](https://quarto.org) manuscript for work on the Sheffield Thyroid Cancer Study
which is the PhD work of [Ovie Edafe (mdp21oe)](https://github.com/mdp21oe).


## GitHub Pages

Because the raw data can not be included in this repository the manuscript can not be rendered via a GitHub Action.

In order to publish the pages you therefore need to do so locally.

1. Complete work on branches.
2. Check the pages render correctly using the _Render_ button in RStudio (shown when you are editing the `index.qmd`
   file).
3. Git commit and push changes to GitHub.
4. Create a Pull Request and merge the changes into the `main` branch.
5. Switch to the `main` branch locally and pull the merged changes down.
6. In the terminal run the following...

``` bash
quarto publish
```

7. Follow the prompts and the site should be rendered and pushed to the `gh-pages` which is setup to be the basis of the
   GitHub hosted pages which can be viewed at
   [ns-rse.github.io/sheffield-thyroid](https://ns-rse.github.io/sheffield-thyroid/). This link can be shared with
   others to see the state of progress.

## Handover 2024-10

### Statistical Analysis

Currently there is code in place that...

1. Cleans and tidies the clinical data.
2. Uses three different methods (Classification And Regression Trees (CART), Predictive Mean Modelling (PMM) and Random
   Forests (RF)) to impute five data sets with each method.
3. Analyse the 15 imputed data sets using the statistical methods of...

  - LASSO Regression
  - Elastic Net Regression
  - Random Forest
  - Gradient Boosting

  ...this is done by splitting each imputed data set into `train` and `test` data, fitting each of the models on the
  training data set and then making predictions in the test cohort to assess how accurate they are.
4. Summaries of these imputed data sets are then available and code in place which aggregates the importance of
   predictors and performance metrics (sensitivity, specificity, positive predictive value, negative predictive value)
   as well as ROC curves in the training and test dataset

### Structure of code

The code is written in various text files in the [Quarto][quarto] literate programming language which interweaves
[Markdown][markdown] with code chunks in the [R][r] programming language.

A lot of the [R][r] code uses the [Tidyverse][tidyverse] and [Tidymodels][tidymodels] frameworks which are opinionated
but coherent approaches to writing R code (they make the code easier to read and reason about compared to base R
commands).

The following sections describe the structure of the work.

#### `index.qmd`

This is the master document. It starts with a [YAML][yaml] header which defines the document structure, authors and how the
code should be compiled (although some settings and configuration options are also in the `_quarto.yaml` file in the
same folder).

It consists of Markdown and code chunks which load packages, describe the data and research. In order to keep the size of
this and other files to a reasonable length not everything is included in this document. Instead there are sections
which call other files, mainly those in the `sections/` directory. Where you see a section as shown below it means the
file listed, the second argument, is included when compiling the document.

``` markdown
{{< include sections/_data_description.qmd >}}
```

Some child documents/files were written and then replaced by newer versions. The code to include these have been left in
place but they have been commented out by placing `<!--` and `-->` around the call to include the document as shown below.


``` markdown
<!-- {{< include sections/_data_description.qmd >}} -->
```

Importantly the first code chunk loads _all_ the libraries used in the analysis (if you introduce a new library then add
it here) as well as setting the current working directory `base_dir` and a number of other directories relative to this
(such as where data is stored in different formats).

#### Cleaning and tidying data

- **File** : [`r/clean.R`](r/clean.R)

This file was written by @mdp21oe with tweaks from @ns-rse and well commented (lines beginning with `##`) to explain
what each section is doing. It saves data to `data/r/df_complete.rds` which can be loaded and used subsequently.

#### Data Description

- **File** : [`sections/_data_description.R`](sections/_data_description.qmd)

This file summarises the data under various sections, it uses a "panel tabset" so that when rendered there are tabs for
each group of characteristics...

- Demographics
- Clinical Characteristics
- Biomarkers
- Ultrasound
- BTA U
- Thyroid Classification
- Cytology

#### Missing Data

- **File** : [`sections/_missing.qmd`][sections/_missing.qmd]

This file summarises the patterns of missing data which are substantial and problematic (more on this below). It
produces tabular and graphical summaries of the missing data, laid out in tabs that group variables logically.


#### Imputation

- **File** : [`sections/_imputation.qmd`](sections/_imputation.qmd)

Because of the large amount of missing data across variables and individuals early attempts at modelling reduced both the
number of observations included and the number of variables in the models. One approach to dealing with missing data is
to impute data, that is to make a prediction of what the missing value should be based on the observed data. This file
undertakes the imputation using the [`mice`](https://amices.org/mice/) package which performs Multivariate Imputation by
Chained Equations. There are many different methods of imputation available and three have been chosen to investigate,
Classification And Regression Trees (CART), Predictive Mean Matching (PMM) and Random Forests (RF). For each method of
imputation five data sets have been imputed. Graphical displays of the distribution of the original and imputed data are
presented which help demonstrate that the imputation method has reasonably estimated the values and not skewed the
distribution.

To reduce the amount of code that needs to be written a function has been introduced (`impute_data()`) which is a
wrapper to imputing the data via the user specified `method` that then produces the tables and plots that summarise the
imputed data.

Each method of imputation saves the five imputed data sets to `data/r/mice_[cart|pmm|rf].rds` for subsequent analysis.

**NB** Personally I am not a great fan of imputation but it is considered a valid and valuable method.

#### Modelling

- **File** : [`sections/_modelling.qmd`](sections/_modelling.qmd)

This section first sources the file [`r/functions.R`](r/functions.R) which defines two functions.

- `model_fitting()` - Fits the a model to a dataset, allows selection of which type of model to fit and which variables
  to include and returns a list of results. If run with `purrr::map()` each element of the list is itself a list of the
  results of having analysed the data set.
- `tidy_imputed_models()` - Extracts data from the a list of imputed analyses (using `purrr:map()`) and extracts tables
  and figures.

The first of these functions (`model_fitting()`) is then used within
[`sections/_modelling.qmd`](sections/_modelling.qmd) in conjunction with [`purrr::map()`]() to analyse each of the
imputed datasets using each of the different modelling approaches. The results are saved to
`data/imputed_[cart|pmm|rf]_model_[lasso|elastic|forest|xgboost].rds` (that is there is one for each imputation method
`cart`, `pmm` or `rf`) and model (`lasso`, `elastic`, `forest` or `xgboost`).

If you want to change the variables included in a model you do so by modifying the call to `model_fitting()`. The
following example reduces the biomarker assays that are included and removes some of the categorical variables.

``` r
## Remove the 'Original' data set from the results of impuation
cart_imputed <- mice_cart$imputed |>
    dplyr::filter(.imp != "Original") |>
    dplyr::mutate(.imp = droplevels(.imp))
## Explicitly state what 'continuous' and which 'categorical' variables are included in the model, these differ from
## the defaults that are set when the function is defined in 'r/functions.R'.
imputed_cart_model_lasso <- cart_imputed |>
    split(cart_imputed$.imp) |>
    purrr::map(\(df) model_fitting(df,
    ## furrr::future_map(\(df) model_fitting(df,
        model_type = "elastic",
        mixture = 1.0,
        imputed = TRUE,
        continuous = c("albumin", "lymphocytes", "monocyte", "size_nodule_mm"),
        categorical = c(
                  "incidental_nodule",
                              "palpable_nodule",
                              "rapid_enlargement",
                              "compressive_symptoms",
                              "hypertension",
                              "vocal_cord_paresis",
                              "graves_disease",
                              "family_history_thyroid_cancer",
                              "bta_u_classification",
                              "solitary_nodule",
                              "cervical_lymphadenopathy",
                              "thy_classification"
                          ),

    ))
tidy_cart_lasso <- tidy_imputed_models(imputed_cart_model_lasso, model_type="elastic")

```

These sections will take a long time to re-run because of the nature of the model fitting and the fact that they are
being run on five data sets. These have an option set to `#| cache: true` which means, in theory, they should only get
re-run when rendering the document if they, or preceding sections
(e.g. [`sections/_imputation.qmd`](sections/_imputation.qmd)) have changed and the results would be impacted.

#### Results

- **Files** :
  - [`sections/_lasso_imputed.qmd`](sections/_lasso_imputed.qmd)
  - [`sections/_elastic_imputed.qmd`](sections/_elastic_imputed.qmd)
  - [`sections/_forest_imputed.qmd`](sections/_forest_imputed.qmd)
  - [`sections/_xgboost_imputed.qmd`](sections/_xgboost_imputed.qmd)

There then follows a series of files which use the results from the [`sections/_modelling.qmd`](sections/_modelling.qmd)
that were saved. The logic being that the results can be saved once and the way in which they are presented can be
modified without having to run the model fitting again.


### Rendering

Currently I've not been able to get the whole document to render correctly. Its a similar error to that encountered with
the multicentre thyroid

### Extracting R code

It can be cumbersome navigating multiple files to run everything in order. A useful strategy to take is to develop code
in the `tmp/scratch.R` file that at the top contains a section which loads all necessary libraries, sets key variables
(as is done by the initial code block of `index.qmd`) and loads all the derived datasets which can then be used. Once
code is working it can be moved to a file under `sections/_<filename>.qmd`. For example when playing around with the
`model_fitting()` function to narrow down subsets you may want to do so in `tmp/scratch.R` and then once you've found
what you want to use update the `sections/_modelling.qmd`.

It is possible to extract all the existing code from the individual `.qmd` files to `.R` files of the same name using
the [`knitr::purl()`](https://bookdown.org/yihui/rmarkdown-cookbook/purl.html) function (it works on Quarto as well as
the described example of RMarkdown). All code must be run which means that child documents that are commented out with
`<!--` and `-->` around them will not be included, you have to un-comment them to have them included. Purl produces one
corresponding `.R` file for each `.qmd` file that is executed (i.e. `sections/_data_description.qmd` produces
`sections/_data_description.R`) rather than producing a single file of all R code.

Further the options at the start of each code chunk includes a directive of whether to explicitly run `purl` on that
code chunk in the form of `#| purl: true`. By default this is `true` for all code chunks and so running...

``` r
knitr::purl("index.qmd")
```

...will result in a `.R` file for every `.qmd`. These could be concatenated using the Linux command `cat` (or the
windows equivalent `type`) although careful attention should be paid to ensure the order of the files


## Challenges and Outstanding work

The biggest challenge with the work is the amount of missing data which is substantial. Part of this is because of bias,
for example certain procedures such as blood assays or biopsies might only be undertaken if there is a strong suspicion
of malignancy.

Because of the amount of missing data and the use of imputation I would be **VERY** cautious about over-interpreting the
results on the sensitivity/specificity/PPV/NPV/AUC ROC statistics.

### Predictor Variables

- Many of the variables are I feel quite likely to be predictors of accuracy as they essential encompass aspects of the
tumor and a clinicians assessment (e.g. `bta_u_classification` and `thy_classification`).
- You may want to tinker with what variables to include in the model in light of how much missing data there is.

### Averaging the imputed models

- For each modelling approach the models fitted to the five imputed datasets of each type need averaging (within each
  imputation approach) and the summary metrics calculating using the averaged model. Not entirely sure how to do this,
  but read that it is done so weighting by the log-likelihood of the model, I wasn't able to work out how to get the
  log-likelihood from each of the modelling approaches, if/when you do add it to the `model_fitting()` function and
  return as an element of the `results`.

### Sensitivity

- Use the `results$fit` to make predictions on the original dataset to get a feel for how effective the model is.

[ovie]: mailto:oedafe1@sheffield.ac.uk
[markdown]: https://www.markdownguide.org/
[quarto]: https://quarto.org
[r]: https://www.r-project.org/
[tidymodels]: https://www.tidymodels.org/
[tidyverse]: https://www.tidyverse.org/
