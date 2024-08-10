## mice imputation trials
df_imp <- mice(df, printFlag = F, seed = 123)
attributes(df_imp)

## complete
library(VIM)
df_trial <- complete(df_imp)
aggr(df_trial)
view(df_trial)
boxplot(df$tsh_value)
boxplot(df_trial$tsh_value)
t.test(df$tsh_value, df_trial$tsh_value)
plot(density(df$tsh_value, na.rm = TRUE), main = "Data wit NA")
lines(density(df_trial$tsh_value, na.rm = TRUE), col = "red", lty = 3)
