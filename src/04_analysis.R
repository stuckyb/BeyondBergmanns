
model_results

# ratio of partial R2 ~ predictors ----
model_results = mutate(model_results, 
                       partial_r2_bio1_bio4 = ifelse(is.na(partial_r2_bio1_bio4), 
                                                     partial_r2_bio1_2yr, 
                                                     partial_r2_bio1_bio4),
                       ratio_temp_precip_r2 = partial_r2_bio1_bio4/partial_r2_bio12_bio15,
                       diff_temp_precip_r2 = partial_r2_bio1_bio4 - partial_r2_bio12_bio15)
					   
hist(model_results$ratio_temp_precip_r2)
hist(model_results$diff_temp_precip_r2)

ggplot(model_results, aes(x = partial_r2_bio1_2yr, fill = taxa)) +
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
paste(predictors, collapse = " + ")

summary(lm(ratio_temp_precip_r2 ~ mass_mean + bio1_mean, data = model_results))
summary(lm(diff_temp_precip_r2 ~ mass_mean + bio1_mean, data = model_results))

mod_1 = lm(diff_temp_precip_r2 ~ mass_mean + bio1_mean + bio1_breadth + bio4_mean + 
             bio4_breadth + bio12_mean + bio12_breadth + bio15_mean + bio15_breadth + 
             habitat_axis1 + habitat_axis2 + trophic_axis1 + trophic_axis2, 
           data = filter(model_results, taxa == "birds"))
summary(step(mod_1))
AIC(lm(diff_temp_precip_r2 ~ bio4_mean + bio15_mean, 
           data = filter(model_results, taxa == "birds")))
summary(lm(diff_temp_precip_r2 ~ bio4_mean + bio15_mean, 
           data = filter(model_results, taxa == "birds")))
# final model
summary(lm(diff_temp_precip_r2 ~ bio4_mean + bio15_mean + habitat_axis2 + trophic_axis1, 
       data = filter(model_results, taxa == "birds")))

mod_2 = lm(diff_temp_precip_r2 ~ mass_mean + bio1_mean + bio1_breadth + bio4_mean + 
             bio4_breadth + bio12_mean + bio12_breadth + bio15_mean + bio15_breadth + 
             habitat_axis1 + habitat_axis2 + trophic_axis1 + trophic_axis2, 
           data = filter(model_results, taxa == "mammals"))
summary(step(mod_2))
# final model
summary(lm(formula = diff_temp_precip_r2 ~ bio1_mean + bio15_mean, 
           data = filter(model_results, taxa == "mammals")))

# elton traits for birds and mammals are generated separately, so they are not comparable.

# coefs ~ predictors ----