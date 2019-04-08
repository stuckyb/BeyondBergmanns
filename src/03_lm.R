# source('src/00_data_funcs.r')
# source('src/01_generate_diagnostics.r')

# separate LMs by sp ----
re_run = F # do you want to re-run all LMs?
if(re_run | (!file.exists(here('data_output/model_results.csv')))){
  model_decisions = read_csv(here('data_raw/metamodeling_data.csv'))
  names(model_decisions)
  model_decisions2 = select(model_decisions, sp = SpName, model = modelName, starts_with('include'))
  lm_by_sp = mutate(model_decisions2,
                    lm_mod = map2(.x = sp, .y = model, .f = lm_1_sp,
                                  data_path = here("data_output/cleaned_bodymass_data/"),
                                  mass_transf = "scale_", 
                                  # log10 or scale_log10 are the other options
                                  verbose = T)) %>% 
    unnest()
  
  lm_by_sp2 = select(lm_by_sp, sp, n_row, starts_with('partial'), 
                     r.squared, adj.r.squared, 
                     r.squared.bio1, adj.r.squared.bio1, 
                     p.value, df, df.residual,
                     AIC, starts_with("(Intercept)"),
                     starts_with("season"),
                     starts_with('bio'))
  
  # combine data
  elton_traits = bind_rows(
    mutate(elton_axes_bird, taxa = 'birds'),
    mutate(elton_axes_mammal, taxa = 'mammals')
  ) %>% 
    rename(sp = Species) %>% 
    mutate(sp = str_replace(sp, " ", "_")) 
  model_results = select(model_decisions, sp = SpName, starts_with('include'),
                         mass_mean = Mean_used, starts_with('bio')) %>% 
    left_join(lm_by_sp2, by = 'sp') %>% 
    left_join(elton_traits, by = 'sp')
  write_csv(model_results, here('data_output/model_results.csv'))
} else {
  model_results = read_csv(here('data_output/model_results.csv'))
}

# hist(model_results$adj.r.squared)

par_r2 = select(model_results, sp, starts_with("partial")) %>% 
  gather('var', 'partial_r2', -sp) %>% 
  mutate(var = recode(var, 
                      'partial_r2_bio1' = 'Temp', 'partial_r2_bio4' = 'Temp S.',
                      'partial_r2_bio12' = 'Precip', 'partial_r2_bio15' = 'Precip S.')) %>% 
  filter(!grepl("partial", var))
var_order = group_by(par_r2, var) %>% summarise(median_r2 = median(partial_r2, na.rm = T)) %>% 
  arrange(desc(median_r2)) %>% 
  pull(var)
par_r2 = mutate(par_r2, var = factor(var, levels = var_order))

# reshape lm coefs
sp_order = arrange(model_results, r.squared)$sp
mod_coefs = gather(model_results, 'var', 'value', starts_with('bio')) %>% 
  separate('var', c('bio_variable', 'var'), sep = "_") %>% 
  select(sp, n_row, r.squared, r.squared.bio1, bio_variable, var, value) %>% 
  spread('var', 'value') %>% 
  mutate(bio_variable = recode(bio_variable, 
                               'bio1' = 'Temp', 'bio4' = 'Temp S.',
                               'bio12' = 'Precip', 'bio15' = 'Precip S.'),
         bio_variable = factor(bio_variable, levels = c('Temp', 'Temp S.', 'Precip', 'Precip S.')),
         sp = factor(sp, levels = sp_order),
         sig = ifelse(p.value < 0.05, "*", ""))

