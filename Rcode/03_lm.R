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
names(lm_by_sp)
head(lm_by_sp)
select(lm_by_sp, sp, multicol_prob, coll_remain, starts_with("partial"), adj.r.squared) %>% View()
select(lm_by_sp, -final_mod, -coef) %>% View()
sum(lm_by_sp$multicol_prob) # 97 sp has potential multicollinearity problem
sum(lm_by_sp$coll_remain) # 5 sp has potential multicollinearity problem after removing bio4

# 36 species with bio4 removed by Narayani
sp_n <- str_replace(c("Aimophila ruficeps", "Aplodontia rufa", "Baeolophus bicolor", "Baiomys taylori", 
                      "Campylorhynchus brunneicapillus", "Chaetodipus formosus", "Chaetodipus penicillatus", 
                      "Chamaea fasciata", "Colinus virginianus", "Corynorhinus townsendii", "Dendroica pinus",
                      "Empidonax occidentalis", "Geomys bursarius", "Melanerpes carolinus", "Melozone fusca",
                      "Microtus longicaudus", "Microtus pennsylvanicus", "Mniotilta varia", "Molothrus ater",
                      "Myodes californicus", "Myotis velifer", "Napaeozapus insignis", "Notiosorex crawfordi",
                      "Oryzomys palustris", "Peromyscus crinitus", "Peromyscus leucopus", "Reithrodontomys fulvescens",
                      "Setophaga ruticilla", "Sigmodon hispidus", "Sorex trowbridgii", "Sorex vagrans", 
                      "Sturnella magna", "Tamias striatus", "Tamiasciurus douglasii", "Thryothorus ludovicianus", 
                      "Tympanuchus phasianellus"), " ", "_")
sp_n2 = filter(lm_by_sp, is.na(partial_r2_bio4))$sp
setdiff(sp_n2, sp_n) # 5
setdiff(sp_n, sp_n2) # 7
length(intersect(sp_n, sp_n2)) # 29 in common

names(lm_by_sp)
lm_by_sp2 = select(lm_by_sp, sp, n_row, multicol_prob, starts_with('partial'), r.squared,
                   adj.r.squared, AIC, starts_with('bio'))
sp_order = arrange(lm_by_sp2, adj.r.squared)$sp
hist(lm_by_sp2$adj.r.squared)

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

library(cowplot)
  
p_a = ggplot(mod_coefs, aes(x = bio_variable, y = sp)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  scale_fill_gradient2() +
  geom_text(aes(label = sig), nudge_y = -0.55) +
  labs(x = 'Climatic predictors', y = '', fill = 'Coefs') +
  theme(legend.position = c(0.44, 0.5),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6.5),
        axis.text.y = element_text(size = 4),
        axis.ticks.y.left = element_line(size = 0.2))

p_b = ggplot(mod_coefs, aes(y = sp, x = adj.r.squared)) +
  geom_segment(aes(yend = sp), xend = 0, color = 'gray') +
  geom_point(color = 'blue', size = 0.8) +
  labs(x = 'Adjust R^2', y = '') +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        plot.margin = margin(l = -0.3, unit = "cm"))
p = plot_grid(p_a, p_b, align = 'h', rel_widths = c(5, 1))
ggsave('figures/lms.pdf', p, width = 8, height = 10)

mod_coefs %>% 
  group_by(bio_variable) %>% 
  summarise(sig_pos = sum(estimate > 0 & p.value < 0.05, na.rm = T),
            sig_neg = sum(estimate < 0 & p.value < 0.05, na.rm = T)) %>% 
  gather("sig", "count", sig_pos, sig_neg) %>% 
  mutate(bio_variable = factor(bio_variable, levels = c('Temp', 'Precip S.', 'Precip', 'Temp S.'))) %>% 
  ggplot(aes(x = bio_variable, y = count, fill = sig)) +
  geom_bar(stat = 'identity', width = 0.5) +
  labs(x = '') +
  coord_flip() +
  scale_fill_manual(values = c('orange', 'blue'))
