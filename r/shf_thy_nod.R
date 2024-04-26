## How to data into R studio, first save file as csv and save into file
## of interest on my computer, then import file into R studio, open a new script
## there will be a code that is provided on the console, copy this to the script
## to read file in R studio, save file into r folder and then ready to go
library(readr)
shf_data <- read_csv("data/csv/sheffield_thyroid_nodule.csv")
## changing names of variables
var_labels <- c(
  study_id = "Study ID", age_at_scan = "Age", ethnicity = "Ethnicity", gender = "Gender",
  ethinicity = "Ethinicity",
  eligibility = "Eligibility", incidental_nodule = "Incidental nodule",
  palpable_nodule = "Palpable nodule", rapid_enlargment = "Rapid enlargement",
  compressive_symtoms = "Compressive symptoms", hypertension = "Hypertension",
  vocal_cord_paresis = "Vocal cord paresis",
  graves_disease = "Graves' disease",
  hashimotos_thyroiditis = "Hashimoto's disease",
  family_history_thyroid_cancer = "Family history of thyroid cancer",
  exposure_radiation = "Exposure to radiation",
  albumin = "Albumin", tsh_value = "TSH value", lymphocytes = "Lymphocytes",
  monocyte = "Monocytes", bta_u_classification = "BTA U", size_nodule_mm = "Nodule size",
  consistency_nodule = "Nodule consistency", cervical_lymphadenopathy = "Cervical Lymphadenopathy",
  repeat_ultrasound = "Repeat ultrasound", repeat_bta_u_classification = "Repeat BTA U",
  fna_done = "FNA done", thy_classification = "Thy classification", repeat_fna_done = "Repeat FNA",
  repeat_thy_classification = "Repeat Thy class", thyroid_surgery = "Thyroid surgery",
  thyroid_histology_diagnosis = "Histology",
  final_pathology = "Final diagnosis", solitary_nodule = "Solitary nodule"
)
## my data set to a data frame so I can compare it to original file (shf_raw)
shf <- tibble(shf_data)
Hmisc::label(shf) <- as.list(var_labels[match(names(shf), names(var_labels))])
## inspect data type for each variable
lapply(shf, typeof)

## converting Yes/No columns into binary 1/0, make a list of all Yes/No Colums

binary_cols <- c(
  "incidental_nodule",
  "palpable_nodule",
  "rapid_enlargment",
  "compressive_symtoms",
  "hypertension",
  "vocal_cord_paresis",
  "graves_disease",
  "hashimotos_thyroiditis",
  "family_history_thyroid_cancer",
  "exposure_radiation",
  "solitary_nodule",
  "cervical_lymphadenopathy",
  "repeat_ultrasound",
  "fna_done",
  "repeat_fna_done"
)

shf <- shf |>
  dplyr::mutate(across(
    all_of(binary_cols),
    ~ dplyr::recode(.x,
      "Yes" = 1,
      "No" = 0
    )
  ))

## converting yes/no variables to variables to logical

shf <- shf |>
  dplyr::mutate(
    incidental_nodule = as.logical(incidental_nodule),
    palpable_nodule = as.logical(palpable_nodule),
    rapid_enlargment = as.logical(rapid_enlargment),
    compressive_symtoms = as.logical(compressive_symtoms),
    hypertension = as.logical(hypertension),
    vocal_cord_paresis = as.logical(vocal_cord_paresis),
    graves_disease = as.logical(graves_disease),
    hashimotos_thyroiditis = as.logical(hashimotos_thyroiditis),
    family_history_thyroid_cancer = as.logical(family_history_thyroid_cancer),
    exposure_radiation = as.logical(exposure_radiation),
    solitary_nodule = as.logical(solitary_nodule),
    cervical_lymphadenopathy = as.logical(cervical_lymphadenopathy),
    repeat_ultrasound = as.logical(repeat_ultrasound),
    fna_done = as.logical(fna_done),
    repeat_fna_done = as.logical(repeat_fna_done)
  )

## converting character variables to factor

shf <- shf |>
  dplyr::mutate(
    gender = as.factor(gender),
    ethnicity = as.factor(ethnicity),
    bta_u_classification = as.factor(bta_u_classification),
    consistency_nodule = as.factor(consistency_nodule),
    repeat_bta_u_classification =
      as.factor(repeat_bta_u_classification),
    thy_classification = as.factor(thy_classification),
    repeat_thy_classification = as.factor(repeat_thy_classification),
    thyroid_surgery = as.factor(thyroid_surgery),
    thyroid_histology_diagnosis = as.factor(thyroid_histology_diagnosis),
    final_pathology = as.factor(final_pathology)
  )

lapply(shf, typeof)

table(shf$gender)

demo_graphics <- select(shf, c(
  "age_at_scan", "gender", "ethnicity"
))

shf |>
  select(c(
    age_at_scan, gender,
    ethnicity,
    final_pathology
  )) |>
  tbl_summary(by = final_pathology)


CreateTableOne(data = shf, vars = c("age_at_scan", "gender", "ethnicity"))

clinical_vars <- select(shf, c(
  "incidental_nodule",
  "palpable_nodule",
  "rapid_enlargment",
  "compressive_symtoms",
  "hypertension",
  "vocal_cord_paresis",
  "graves_disease",
  "hashimotos_thyroiditis",
  "family_history_thyroid_cancer",
  "exposure_radiation"
))
CreateTableOne(data = clinical_vars, vars = c(
  "incidental_nodule",
  "palpable_nodule",
  "rapid_enlargment",
  "compressive_symtoms",
  "hypertension",
  "vocal_cord_paresis",
  "graves_disease",
  "hashimotos_thyroiditis",
  "family_history_thyroid_cancer",
  "exposure_radiation"
))

shf |>
  select(c(
    incidental_nodule,
    palpable_nodule,
    rapid_enlargment,
    compressive_symtoms,
    hypertension,
    vocal_cord_paresis,
    graves_disease,
    hashimotos_thyroiditis,
    family_history_thyroid_cancer,
    exposure_radiation, final_pathology
  )) |>
  tbl_summary(by = final_pathology) |>
  add_p()
## need to check how the above deals with missing data

## biochemical variables summary

biochem_vars <- select(shf, c(
  "albumin",
  "tsh_value",
  "lymphocytes",
  "monocyte"
))

shf |>
  select(c(
    albumin,
    tsh_value,
    lymphocytes,
    monocyte, final_pathology
  )) |>
  tbl_summary(by = final_pathology) |>
  add_p()

CreateTableOne(data = biochem_vars, vars = c(
  "albumin",
  "tsh_value",
  "lymphocytes",
  "monocyte"
))

## imaging and thyc

ultrasound_char <- select(shf, c(
  "size_nodule_mm",
  "solitary_nodule",
  "bta_u_classification",
  "consistency_nodule",
  "cervical_lymphadenopathy"
))

shf |>
  select(c(
    size_nodule_mm, solitary_nodule,
    bta_u_classification, consistency_nodule,
    cervical_lymphadenopathy,
    final_pathology
  )) |>
  tbl_summary(by = final_pathology) |>
  add_p()

CreateTableOne(data = ultrasound_char, vars = c(
  "size_nodule_mm",
  "solitary_nodule",
  "bta_u_classification",
  "consistency_nodule",
  "cervical_lymphadenopathy"
))

cytology_char <- select(shf, c(
  "fna_done",
  "thy_classification",
  "repeat_fna_done",
  "thyroid_surgery",
  "final_pathology"
))

CreateTableOne(data = cytology_char, vars = c(
  "fna_done",
  "thy_classification",
  "repeat_fna_done",
  "thyroid_surgery",
  "final_pathology"
))

shf |>
  select(c(
    thy_classification,
    final_pathology
  )) |>
  tbl_summary(by = final_pathology)

## create dir where file is saved
