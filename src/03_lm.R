# source('src/00_data_funcs.r')
# source('src/01_generate_diagnostics.r')

model_decisions = read_csv(here('data_raw/metamodeling_data.csv'))
names(model_decisions)
model_decisions2 = select(model_decisions, sp = SpName, model = modelName, starts_with('include'))
lm_by_sp = mutate(model_decisions2,
                       lm_mod = map2(.x = sp, .y = model, .f = lm_1_sp,
                                     data_path = here("data_output/cleaned_bodymass_data/"),
                                     mass_transf = "log10", 
                                     # scale_ or scale_log10 are the other options
                                     verbose = T)) %>% 
  unnest()
select(lm_by_sp, -final_mod, -coef) %>% View()                      

lm_by_sp2 = select(lm_by_sp, sp, n_row, starts_with('partial'), r.squared,
                   adj.r.squared, AIC, starts_with('bio'))
sp_order = arrange(lm_by_sp2, adj.r.squared)$sp
hist(lm_by_sp2$adj.r.squared)

par_r2 = select(lm_by_sp2, sp, starts_with("partial")) %>% 
  gather('var', 'partial_r2', -sp) %>% 
  mutate(var = recode(var, 
                      'partial_r2_bio1' = 'Temp', 'partial_r2_bio4' = 'Temp S.',
                      'partial_r2_bio12' = 'Precip', 'partial_r2_bio15' = 'Precip S.'))
var_order = group_by(par_r2, var) %>% summarise(median_r2 = median(partial_r2, na.rm = T)) %>% 
  arrange(desc(median_r2)) %>% 
  pull(var)
par_r2 = mutate(par_r2, var = factor(var, levels = var_order))
  

# reshape lm coefs
mod_coefs = gather(lm_by_sp2, 'var', 'value', starts_with('bio')) %>% 
  separate('var', c('bio_variable', 'var'), sep = "_") %>% 
  select(sp, n_row, adj.r.squared, bio_variable, var, value) %>% 
  spread('var', 'value') %>% 
  mutate(bio_variable = recode(bio_variable, 
                               'bio1' = 'Temp', 'bio4' = 'Temp S.',
                               'bio12' = 'Precip', 'bio15' = 'Precip S.'),
         bio_variable = factor(bio_variable, levels = c('Temp', 'Temp S.', 'Precip', 'Precip S.')),
         sp = factor(sp, levels = sp_order),
         sig = ifelse(p.value < 0.05, "*", "")) 
