names(d)
d
length(unique(d$sp)) # 175 sp
# for each sp, conduct a lm

lm_by_sp = tibble(sp = unique(d$sp))
lm_by_sp = mutate(lm_by_sp,
                  lm_mod = map(sp, lm_1_sp,
                               fm = "mass_g_log10 ~  bio1 + bio12 + bio4 + bio15 + season",
                               df = d)) %>% 
  unnest()
select(lm_by_sp, species, vif_mean, coll_prob) %>% View()
sum(lm_by_sp$coll_prob) # 45 sp has potential multicollinearity problem
