## How to data into R studio, first save file as csv and save into file
## of interest on my computer, then import file into R studio, open a new script
## there will be a code that is provided on the console, copy this to the script
## to read file in R studio, save file into r folder and then ready to go

library(readr)
## ns-rse (2024-04-29) - The following line has been tweaked to work with the environment variable csv_dir which is set
## in ../index.qmd  as this is the file that runs all of the analysis and produces the report. Further the raw data is
## read to `raw_data`. A copy is made further down to `df` and this is cleaned. As the data only pertains to one
## location there is no benefit to having it as `shd`
raw_data <- readr::read_csv(paste(csv_dir, "sheffield_thyroid_nodule.csv", sep = "/"))

## Set labels for each variable
## ns-rse (2024-04-29) : Try and keep each expression in long lists like this on a separate line, where possible having
## them sorted alphabetically really helps make your life easier when you want to find and change something.
var_labels <- c(
  age_at_scan = "Age",
  albumin = "Albumin",
  bta_u_classification = "BTA U",
  cervical_lymphadenopathy = "Cervical Lymphadenopathy",
  compressive_symptoms = "Compressive symptoms",
  consistency_nodule = "Nodule consistency",
  eligibility = "Eligibility",
  ethinicity = "Ethinicity",
  ethnicity = "Ethnicity",
  exposure_radiation = "Exposure to radiation",
  family_history_thyroid_cancer = "Family history of thyroid cancer",
  final_pathology = "Final diagnosis",
  fna_done = "FNA done",
  gender = "Gender",
  graves_disease = "Graves' disease",
  hashimotos_thyroiditis = "Hashimoto's disease",
  hypertension = "Hypertension",
  incidental_nodule = "Incidental nodule",
  lymphocytes = "Lymphocytes",
  monocyte = "Monocytes",
  palpable_nodule = "Palpable nodule",
  rapid_enlargement = "Rapid enlargement",
  repeat_bta_u_classification = "Repeat BTA U",
  repeat_fna_done = "Repeat FNA",
  repeat_thy_classification = "Repeat Thy class",
  repeat_ultrasound = "Repeat ultrasound",
  size_nodule_mm = "Nodule size (mm)",
  solitary_nodule = "Solitary nodule",
  study_id = "Study ID",
  thy_classification = "Thy classification",
  thyroid_histology_diagnosis = "Histology",
  thyroid_surgery = "Thyroid surgery",
  tsh_value = "TSH value",
  vocal_cord_paresis = "Vocal cord paresis"
  )

## @ns-rse 2024-07-18 : I noticed there is a single instance of graves_disease that is "no" rather than "No", we correct
## that now
raw_data <- raw_data|>
    dplyr::mutate(graves_disease = case_when(graves_disease == "no" ~ "No",
        .default = graves_disease
        ))

## @ns-rse 2024-08-13 : Convert is.na(thy_classification) to `No biopsy` and tidy up other variables too.
##                      Based on lines 107-119 of index.qmd from commit :
##                      https://github.com/ns-rse/sheffield-thyroid/commit/9d586f8e701a015e27c693a8c5788dbd74e1ac2d
raw_data <- raw_data|>
    dplyr::mutate(thy_classification = dplyr::case_when(is.na(thy_classification) ~ "No biopsy",
                                                        .default = thy_classification),
                  exposure_radiation = dplyr::case_when(is.na(exposure_radiation) ~ "No",
                                                        .default = exposure_radiation),
                  repeat_fna_done = dplyr::case_when(is.na(repeat_fna_done) ~ "No",
                                                     .default = repeat_fna_done),
                  repeat_ultrasound = dplyr::case_when(is.na(repeat_ultrasound) ~ "No",
                                                       .default = repeat_ultrasound),
                  repeat_bta_u_classification = dplyr::case_when(is.na(repeat_bta_u_classification) ~ "No ultrasound",
                                                                 .default = repeat_bta_u_classification),
                  repeat_thy_classification = dplyr::case_when(is.na(repeat_thy_classification) ~ "No biopsy",
                                                               .default = repeat_thy_classification),
                  thyroid_histology_diagnosis = dplyr::case_when(is.na(thyroid_histology_diagnosis) ~ "No surgery",
                                                                 .default = thyroid_histology_diagnosis))


## Make a copy of the raw data so that comparisons can be made between it and the cleaned dataset. This allows checking
## that no mistakes have been made
df <- tibble(raw_data)
Hmisc::label(df) <- as.list(var_labels[match(names(df), names(var_labels))])


## Convert character variables to factors, this now includes binary variables that are 'No'/'Yes' which will now be
## treated as factors /nominal variables and the required dummies generated when we use
## recipes::step_dummy(recipes::all_nominal_predictors()) to encode them to dummy variables.
df <- df |>
    dplyr::mutate_if(is.character, as.factor)

## ns-rse (2024-04-30) : A Check should be made for duplicated data and any such duplicates removed, even if you think
## there shouldn't be or aren't any duplicates in a dataset it is a simple step to make such a check and remove any that
## are found.
df <- unique(df)

## Save a copy of the clean data (df) for loading and using in analysis
saveRDS(df, file = paste(r_dir, "clean.rds", sep = "/"))

## @ns-rse 2024-08-13 : Moving the select and filter to produce `df_complete` to here and saving as an object for
## loading.
## @ns-rse 2024-06-18:
## We want this table to match the model, therefore rather than repeat ourselves we move the subsetting to give us
## df_complete to here and assign to df_complete. This is then passed into gtsummary::tbl_summary() and is available for
## the next code chunk.
##
## This is the point at which you should subset your data for those who have data available for the variables of
## interest. The variables here should include the outcome `final_pathology` and the predictor variables that are set in
## the code chunk `recipe`. May want to move this to another earlier in the processing so that the number of rows can be
## counted and reported.
df_complete <- df |>
  dplyr::select(
    age_at_scan,
    gender,
    ethnicity,
    incidental_nodule,
    palpable_nodule,
    rapid_enlargement,
    compressive_symptoms,
    hypertension,
    vocal_cord_paresis,
    graves_disease,
    hashimotos_thyroiditis,
    family_history_thyroid_cancer,
    exposure_radiation,
    albumin,
    tsh_value,
    lymphocytes,
    monocyte,
    bta_u_classification,
    solitary_nodule,
    size_nodule_mm,
    cervical_lymphadenopathy,
    thy_classification,
    final_pathology) |>
## @ns-rse 2024-06-14 :
## I would consider removing this dplyr::filter() and instead using the recipes::step_filter_missing() as is now in
## place below
## dplyr::filter(if_any(everything(), is.na))
## Instead I think it might be useful to remove individuals who do not have a value for final_pathology here
  dplyr::filter(!is.na(final_pathology))

## Save a copy of the clean data (df) for loading and using in analysis
saveRDS(df_complete, file = paste(r_dir, "df_complete.rds", sep = "/"))


## OBSOLETE CODE
## The following has been replaced with more succinct versions

## ns-rse (2024-04-29) : Duplicates the above, no point going to 0/1 then FALSE/TRUE just go direct
## Convert binary (yes/no) variables to logical
## df <- df |>
##   dplyr::mutate(
##     incidental_nodule = as.logical(incidental_nodule),
##     palpable_nodule = as.logical(palpable_nodule),
##     rapid_enlargement = as.logical(rapid_enlargement),
##     compressive_symptoms = as.logical(compressive_symptoms),
##     hypertension = as.logical(hypertension),
##     vocal_cord_paresis = as.logical(vocal_cord_paresis),
##     graves_disease = as.logical(graves_disease),
##     hashimotos_thyroiditis = as.logical(hashimotos_thyroiditis),
##     family_history_thyroid_cancer = as.logical(family_history_thyroid_cancer),
##     exposure_radiation = as.logical(exposure_radiation),
##     solitary_nodule = as.logical(solitary_nodule),
##     cervical_lymphadenopathy = as.logical(cervical_lymphadenopathy),
##     repeat_ultrasound = as.logical(repeat_ultrasound),
##     fna_done = as.logical(fna_done),
##     repeat_fna_done = as.logical(repeat_fna_done)
##   )



## Convert character variables to factor
## df <- df |>
##   dplyr::mutate(
##     gender = as.factor(gender),
##     ethnicity = as.factor(ethnicity),
##     bta_u_classification = as.factor(bta_u_classification),
##     consistency_nodule = as.factor(consistency_nodule),
##     repeat_bta_u_classification = as.factor(repeat_bta_u_classification),
##     thy_classification = as.factor(thy_classification),
##     repeat_thy_classification = as.factor(repeat_thy_classification),
##     thyroid_surgery = as.factor(thyroid_surgery),
##     thyroid_histology_diagnosis = as.factor(thyroid_histology_diagnosis),
##     final_pathology = as.factor(final_pathology)
##   )

## REDUNDANT CODE
##
## The following code tabulates variables for checking, its not essential to the cleaning of the data, most tables will
## be summarised and presented in the manuscript.
## table(df$gender)

## lapply(df, typeof)

## demo_graphics <- select(df, c(
##   "age_at_scan", "gender", "ethnicity"
## ))

## df |>
##   select(c(
##     age_at_scan, gender,
##     ethnicity,
##     final_pathology
##   )) |>
##   tbl_summary(by = final_pathology)


## CreateTableOne(data = df, vars = c("age_at_scan", "gender", "ethnicity"))

## clinical_vars <- select(df, c(
##   "incidental_nodule",
##   "palpable_nodule",
##   "rapid_enlargement",
##   "compressive_symptoms",
##   "hypertension",
##   "vocal_cord_paresis",
##   "graves_disease",
##   "hashimotos_thyroiditis",
##   "family_history_thyroid_cancer",
##   "exposure_radiation"
## ))
## CreateTableOne(data = clinical_vars, vars = c(
##   "incidental_nodule",
##   "palpable_nodule",
##   "rapid_enlargement",
##   "compressive_symptoms",
##   "hypertension",
##   "vocal_cord_paresis",
##   "graves_disease",
##   "hashimotos_thyroiditis",
##   "family_history_thyroid_cancer",
##   "exposure_radiation"
## ))

## df |>
##   select(c(
##     incidental_nodule,
##     palpable_nodule,
##     rapid_enlargement,
##     compressive_symptoms,
##     hypertension,
##     vocal_cord_paresis,
##     graves_disease,
##     hashimotos_thyroiditis,
##     family_history_thyroid_cancer,
##     exposure_radiation, final_pathology
##   )) |>
##   tbl_summary(by = final_pathology) |>
##   add_p()
## ## need to check how the above deals with missing data

## ## biochemical variables summary

## biochem_vars <- select(df, c(
##   "albumin",
##   "tsh_value",
##   "lymphocytes",
##   "monocyte"
## ))

## df |>
##   select(c(
##     albumin,
##     tsh_value,
##     lymphocytes,
##     monocyte, final_pathology
##   )) |>
##   tbl_summary(by = final_pathology) |>
##   add_p()

## CreateTableOne(data = biochem_vars, vars = c(
##   "albumin",
##   "tsh_value",
##   "lymphocytes",
##   "monocyte"
## ))

## ## imaging and thyc

## ultrasound_char <- select(df, c(
##   "size_nodule_mm",
##   "solitary_nodule",
##   "bta_u_classification",
##   "consistency_nodule",
##   "cervical_lymphadenopathy"
## ))

## df |>
##   select(c(
##     size_nodule_mm, solitary_nodule,
##     bta_u_classification, consistency_nodule,
##     cervical_lymphadenopathy,
##     final_pathology
##   )) |>
##   tbl_summary(by = final_pathology) |>
##   add_p()

## CreateTableOne(data = ultrasound_char, vars = c(
##   "size_nodule_mm",
##   "solitary_nodule",
##   "bta_u_classification",
##   "consistency_nodule",
##   "cervical_lymphadenopathy"
## ))

## cytology_char <- select(df, c(
##   "fna_done",
##   "thy_classification",
##   "repeat_fna_done",
##   "thyroid_surgery",
##   "final_pathology"
## ))

## CreateTableOne(data = cytology_char, vars = c(
##   "fna_done",
##   "thy_classification",
##   "repeat_fna_done",
##   "thyroid_surgery",
##   "final_pathology"
## ))

## df |>
##   select(c(
##     thy_classification,
##     final_pathology
##   )) |>
##   tbl_summary(by = final_pathology)

## ## create dir where file is saved
