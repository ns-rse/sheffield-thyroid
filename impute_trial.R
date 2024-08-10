## mice imputation trials
df_imp <- mice(df, printFlag = F, seed = 123)
attributes(df_imp)

## complete
library(VIM)
df_imp <- complete(df_imp)
aggr(df_imp)
view(df_imp)
boxplot(df$tsh_value)
boxplot(df_imp$tsh_value)
t.test(df$tsh_value, df_imp$tsh_value)
       