source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/plot_res.R')
df <- df %>% 
  select(-patient_id, -sample, -sample_class, -cluster, -RMSE, -cna_caller, -mut_caller) %>% 
  mutate(type = factor(type, levels = c('SBS', 'ID')),
         tool = factor(tool, levels = c('viber', 'pyclonevi', 'mobster_univariate')),
         sig_tool = factor(sig_tool, levels = c('SigProfiler', 'BASCULE')))


# Running the linear regression
model <- lm(CosineSimilarity ~ ., data = df)

# Viewing the results
summary(model)



# Tidy the model results
model_tidied <- tidy(model, conf.int = TRUE) |> 
  filter(term != "(Intercept)") # Usually excluded to keep the scale readable

# Create the plot
ggplot(model_tidied, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  labs(title = "Effect of Predictors on CosineSimilarity",
       subtitle = "Dots represent the Estimate; bars represent 95% Confidence Intervals",
       x = "Estimate (Change in CosineSimilarity)",
       y = "Predictor Variable") +
  theme_minimal()


par(mfrow = c(2, 2))
plot(model)

