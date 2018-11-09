
model_results

# ratio of partial R2 ~ predictors ----
model_results = mutate(model_results, 
                       partial_r2_bio1_bio4 = ifelse(is.na(partial_r2_bio1_bio4), 
                                                     partial_r2_bio1, 
                                                     partial_r2_bio1_bio4),
                       ratio_temp_precip_r2 = partial_r2_bio1_bio4/partial_r2_bio12_bio15,
                       diff_temp_precip_r2 = partial_r2_bio1_bio4 - partial_r2_bio12_bio15)
hist(model_results$ratio_temp_precip_r2)
hist(model_results$diff_temp_precip_r2)
ggplot(model_results, aes(x = partial_r2_bio1, fill = taxa)) +
  geom_histogram() +
  facet_wrap(~taxa, ncol = 1)

model_results = mutate(model_results, mass_mean = log10(mass_mean))

predictors = c("mass_mean", 
               "bio1_mean", "bio1_breadth", 
               "bio4_mean", "bio4_breadth", 
               "bio12_mean", "bio12_breadth",
               "bio15_mean", "bio15_breadth",
               "habitat_axis1", "habitat_axis2", 
               "trophic_axis1", "trophic_axis2")


summary(lm(ratio_temp_precip_r2 ~ mass_mean + bio1_mean, data = model_results))
summary(lm(diff_temp_precip_r2 ~ mass_mean + bio1_mean, data = model_results))

mod_1 = lm(diff_temp_precip_r2 ~ mass_mean + bio1_mean + bio1_breadth + bio12_mean + bio12_breadth, 
           data = filter(model_results, taxa == "birds"))
summary(mod_1)

# elton traits for birds and mammals are generated separately, so they are not comparable.

# coefs ~ predictors ----