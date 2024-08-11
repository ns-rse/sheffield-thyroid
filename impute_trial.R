## mice imputation trials

##df_imp <- mice(df,
#               m = 5,
#               method = Null
#               )

##method(mice)
#df_imp <- mice(df, printFlag = F, seed = 123)
#attributes(df_imp)

## complete
#library(VIM)
#df_trial <- complete(df_imp)
#aggr(df_trial)
#view(df_trial)
#boxplot(df$tsh_value)
#boxplot(df_trial$tsh_value)
#t.test(df$tsh_value, df_trial$tsh_value)
#plot(density(df$tsh_value, na.rm = TRUE), main = "Data wit NA")
#lines(density(df_trial$tsh_value, na.rm = TRUE), col = "red", lty = 3)

## another take on imputation with first visualing missingness
## and selecting specific variable to include and methods for each variable

library(mice)
library(tidyverse)
library(naniar)
library(VIM)
summary(df)

## first, visualise missingness

## Visualize the overall missingness in the dataset
vis_miss(df)

## Visualize missing data by variable, 11.1% missing data now, after recoding some
## of the NA values, nodule consistency was infrequently
## reported by the radiologist, too much of this missing to be imputated i feel > 50%
gg_miss_var(df)

## Visualize missing data by case (rows), nodule consistency was infrequently
## reported by the radiologist, too much of this missing to be imputated i feel
gg_miss_case(df)

## Visualize the pattern of missing data using an upset plot, interesting plot
## looking at combination of missingness, makes sense the bloods tests frequently
## "missing" but better termed not done. Missing at Random here.
gg_miss_upset(df)

## below plot not very useful as data set too big
aggr(df, col = c('navyblue','red'),
     numbers = TRUE, sortVars = TRUE, labels = names(df),
     cex.axis = .8, gap = 3, ylab = c("Missing data","Pattern"))

## code to exclude study_id and eligibility before viewing missingness

exclude_vars <- c("study_id", "eligibility")
missing_data <- data.frame(
  variable = setdiff(names(df), exclude_vars),
  missing = colSums(is.na(df))[setdiff(names(df), exclude_vars)]
)


### this plot visualise missingness by variable
ggplot(missing_data, aes(x = reorder(variable, -missing), y = missing)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Missing Data by Variable", x = "Variable", y = "Number of Missing Values") +
  coord_flip()

# select which variables to include in the imputation, specifically removed
## study ID, final pathology, and other variables not used for prediction
# Create the predictor matrix with the included variables
pred_matrix <- mice::quickpred(df, 
                               include =  c("age_at_scan",
                                           "gender",
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
                                           "albumin",
                                           "tsh_value",
                                           "lymphocytes",
                                           "monocyte",
                                           "bta_u_classification",
                                           "solitary_nodule",
                                           "size_nodule_mm",
                                           "cervical_lymphadenopathy",
                                           "thy_classification"))

## Manually list variables to exclude from the predictor matrix and imputation

vars_to_exclude <- c("study_id", 
                     "eligibility", 
                     "repeat_ultrasound",
                     "fna_done",
                     "repeat_bta_u_classification",
                     "repeat_thy_classification",
                     "final_pathology")

pred_matrix[vars_to_exclude, ] <- 0
pred_matrix[, vars_to_exclude] <- 0  

## View the predictor matrix, for some reason the print still has variable excluded, unsure if this is still being used
## for imputation
print(pred_matrix)

## specify methods to use for imputation for each variable

## Define methods based on the variable types
method <- make.method(df)
## I noted R can automatically select the appropraite method of imputation, but here I have
## specified the method for each variable manually
method["age_at_scan"] <- "pmm"
method["gender"] <- "logreg"
method["ethnicity"] <- "polyreg"
method["incidental_nodule"] <- "logreg"
method["palpable_nodule"] <- "logreg"
method["rapid_enlargement"] <- "logreg"
method["compressive_symptoms"] <- "logreg"
method["hypertension"] <- "logreg"
method["vocal_cord_paresis"] <- "logreg"
method["graves_disease"] <- "logreg"
method["hashimotos_thyroiditis"] <- "logreg"
method["family_history_thyroid_cancer"] <- "logreg"
method["exposure_radiation"] <- "logreg"
method["albumin"] <- "pmm"
method["tsh_value"] <- "pmm"
method["lymphocytes"] <- "pmm"
method["monocyte"] <- "pmm"
method["bta_u_classification"] <- "polyreg"
method["solitary_nodule"] <- "logreg"
method["size_nodule_mm"] <- "pmm"
method["cervical_lymphadenopathy"] <- "logreg"
method["thy_classification"] <- "polyreg"

## Perform imputation on the pred matrix I created above
imputed_data <- mice(df, 
                     method = method, 
                     predictorMatrix = pred_matrix, 
                     m = 5,
                     maxit = 10,
                     seed = 123)
summary(imputed_data)

## this gives the complete imputed data
complete_data <- complete(imputed_data)

summary(complete_data)

## to compare missing values in original data set to imputed data, i will
## first identify the various imputations
imputed_data_1 <- complete(imputed_data, action = 1)
imputed_data_2 <- complete(imputed_data, action = 2)
imputed_data_3 <- complete(imputed_data, action = 3)
imputed_data_4 <- complete(imputed_data, action = 4)
imputed_data_5 <- complete(imputed_data, action = 5)

## now I will try to visualise the imputations in comparison with the original data
## using tsh  value

df_combined <- data.frame(
  original = df$tsh_value,
  imputed_1 = imputed_data_1$tsh_value,
  imputed_2 = imputed_data_2$tsh_value,
  imputed_3 = imputed_data_3$tsh_value,
  imputed_4 = imputed_data_4$tsh_value,
  imputed_5 = imputed_data_5$tsh_value
)

## visualising data frame for imputed data TSH value versus original data &
## using t test to evaluate if there is a significant difference between imputated and 
## original data set

df_combined <- data.frame(
  value = c(df$tsh_value, imputed_data_1$tsh_value),
  dataset = rep(c("Original", "Imputed"), each = nrow(df))
)
ggplot(df_combined, aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot() +
  labs(title = "Boxplot of TSH Values: Original vs Imputed Data",
       x = "Dataset",
       y = "TSH Value") +
  theme_minimal()


mean(df$tsh_value, na.rm = TRUE)
mean(imputed_data_1$tsh_value)
t.test(df$tsh_value, imputed_data_1$tsh_value)

df_combined <- data.frame(
  value = c(df$tsh_value, imputed_data_2$tsh_value),
  dataset = rep(c("Original", "Imputed"), each = nrow(df))
)
ggplot(df_combined, aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot() +
  labs(title = "Boxplot of TSH Values: Original vs Imputed Data",
       x = "Dataset",
       y = "TSH Value") +
  theme_minimal()
mean(df$tsh_value, na.rm = TRUE)
mean(imputed_data_2$tsh_value)
t.test(df$tsh_value, imputed_data_1$tsh_value)

df_combined <- data.frame(
  value = c(df$tsh_value, imputed_data_3$tsh_value),
  dataset = rep(c("Original", "Imputed"), each = nrow(df))
)
ggplot(df_combined, aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot() +
  labs(title = "Boxplot of TSH Values: Original vs Imputed Data",
       x = "Dataset",
       y = "TSH Value") +
  theme_minimal()
mean(df$tsh_value, na.rm = TRUE)
mean(imputed_data_3$tsh_value)
t.test(df$tsh_value, imputed_data_3$tsh_value)

df_combined <- data.frame(
  value = c(df$tsh_value, imputed_data_4$tsh_value),
  dataset = rep(c("Original", "Imputed"), each = nrow(df))
)
ggplot(df_combined, aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot() +
  labs(title = "Boxplot of TSH Values: Original vs Imputed Data",
       x = "Dataset",
       y = "TSH Value") +
  theme_minimal()
mean(df$tsh_value, na.rm = TRUE)
mean(imputed_data_4$tsh_value)
t.test(df$tsh_value, imputed_data_4$tsh_value)

df_combined <- data.frame(
  value = c(df$tsh_value, imputed_data_5$tsh_value),
  dataset = rep(c("Original", "Imputed"), each = nrow(df))
)
ggplot(df_combined, aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot() +
  labs(title = "Boxplot of TSH Values: Original vs Imputed Data",
       x = "Dataset",
       y = "TSH Value") +
  theme_minimal()
mean(df$tsh_value, na.rm = TRUE)
mean(imputed_data_5$tsh_value)
t.test(df$tsh_value, imputed_data_5$tsh_value)


