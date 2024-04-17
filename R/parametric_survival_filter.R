







#####The exponential 
library(epiR)
library(ggsurvfit)
s1 <- survfit(Surv(days, vital_status) ~ 1, data = survival)

survfit2(Surv(days, vital_status) ~ 1, data = survival) %>%  #survival plot
  ggsurvfit(color = "#73a9c9", size = 1.2) +
  labs(
    title = "TCGA Kaplan-Meier Survival Function",
    x = "Days",
    y = "Overall survival probability"
  )


epiR::epi.insthaz(s1) %>%  #Hazard function plot
  ggplot(aes(x = time, y = hest), color = 'red') + 
  geom_smooth(color = "red", method = "loess", formula = "y ~ x") +
  theme_light() +
  labs(title = "Kaplan-Meier Hazard Function Estimate", 
       x = "Time", y = "Instantaneous Hazard")

# par_fits <- tibble(
#   dist_param = c("exp", "weibull", "gompertz", "gamma", "lognormal", "llogis",
#                  "gengamma"),
#   dist_name = c("Exponential", "Weibull (AFT)", "Gompertz", "Gamma",
#                 "Lognormal", "Log-logistic", "Generalized gamma")


par_fits <- tibble(
  dist_param = c("exp", "weibull", "gamma", "lognormal", "llogis","gengamma"),
  dist_name = c("Exponential", "Weibull (AFT)", "Gamma", "Lognormal", "Log-logistic", "Generalized gamma")
) %>%
  mutate(
    fit = map(dist_param, ~flexsurvreg(Surv(days, vital_status) ~ 1, data = survival, dist = .x)),
    fit_smry = map(fit, ~summary(.x, type = "hazard", ci = TRUE, tidy = TRUE)),
    AIC = map_dbl(fit, ~.x$AIC)
  )


plot_par_fits <- par_fits %>%
  select(-c(dist_param, fit)) %>%
  unnest(fit_smry)


ggplot() +
  geom_line(data = plot_par_fits, aes(x = time, y = est, color = dist_name)) +
  geom_smooth(data = epiR::epi.insthaz(s1), method = "loess", formula = "y ~ x", aes(x = time, y = hest), color = 'red') +
  theme_light() +
  labs(title = "Parametric Distribution Fits to TCGA Data.")

par_fits %>%
  arrange(AIC) %>%
  select(dist_name, AIC)




#Create the synthlethal parameter for every pair


exp.mod.aft <- survreg(Surv(days, vital_status) ~ gender + years_to_birth, data = survival, dist = 'exponential')
summary(exp.mod.aft)

SL_pairs_binarize_expression_path <-
  file.path(datapath, 'SL_pairs_binarize_expression.csv')
SL_pairs <- read_csv(SL_pairs_binarize_expression_path)

for (i in  1:nrow(SL_pairs)){
  if (SL_pairs$gene1[i] %in% expression$gene && SL_pairs$gene2[i] %in% expression$gene){
    
    
    survival <- survival |> mutate(SL_index = ) 
    survreg(Surv(days, vital_status) ~ gender + years_to_birth + SL_index, data = survival, dist = 'exponential')$coefficients["SL_index"]
    
  }
}
