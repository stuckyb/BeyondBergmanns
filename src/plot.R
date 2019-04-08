# Fig 1 ----
# mod_ceof2 from 03_lm.R
sp_taxa = data.frame(sp = factor(sp_order, levels = sp_order)) %>% 
  left_join(select(model_results, sp, taxa) %>% 
              mutate(sp = factor(sp, levels = sp_order))) %>% 
  mutate(co = ifelse(taxa == "birds", "olivedrab", "black"))
p_a_birds = ggplot(filter(mod_coefs, sp %in% filter(sp_taxa, taxa == "birds")$sp), 
             aes(x = bio_variable, y = sp)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  scale_fill_gradient2() +
  geom_text(aes(label = sig), nudge_y = -0.3) +
  labs(x = 'Climatic predictors', y = '', fill = 'Coefs') +
  theme(legend.position = c(0.44, 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.y.left = element_line(),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))
p_b_birds = ggplot(filter(mod_coefs, sp %in% filter(sp_taxa, taxa == "birds")$sp), 
             aes(y = sp, x = r.squared)) +
  geom_segment(aes(yend = sp), xend = 0, color = 'gray') +
  geom_segment(aes(yend = sp, x = r.squared.bio1, y = sp), xend = 0, 
               color = 'mediumvioletred', inherit.aes = F) +
  geom_point(color = 'blue', size = 0.8) +
  labs(x = expression(paste("", R^{2})), 
       y = '') +
  xlim(c(0, 0.5)) +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        plot.margin = margin(l = -0.3, unit = "cm"),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))
p_birds = plot_grid(p_a_birds, p_b_birds, align = 'h', rel_widths = c(5, 1))
ggsave(here('figures/fig1_lms_birds.png'), p_birds, width = 9, height = 11)
ggsave(here('figures/fig1_lms_birds.pdf'), p_birds, width = 9, height = 11)

p_a_mammals = ggplot(filter(mod_coefs, sp %in% filter(sp_taxa, taxa == "mammals")$sp), 
                   aes(x = bio_variable, y = sp)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  scale_fill_gradient2() +
  geom_text(aes(label = sig), nudge_y = -0.3) +
  labs(x = 'Climatic predictors', y = '', fill = 'Coefs') +
  theme(legend.position = c(0.44, 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.y.left = element_line(),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))
p_b_mammals = ggplot(filter(mod_coefs, sp %in% filter(sp_taxa, taxa == "mammals")$sp), 
                   aes(y = sp, x = r.squared)) +
  geom_segment(aes(yend = sp), xend = 0, color = 'gray') +
  geom_segment(aes(yend = sp, x = r.squared.bio1, y = sp), xend = 0, 
               color = 'mediumvioletred', inherit.aes = F) +
  geom_point(color = 'blue', size = 0.8) +
  labs(x = expression(paste("", R^{2})), 
       y = '') +
  xlim(c(0, 0.5)) +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        plot.margin = margin(l = -0.3, unit = "cm"),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))
p_mammals = plot_grid(p_a_mammals, p_b_mammals, align = 'h', rel_widths = c(5, 1))
ggsave(here('figures/fig1_lms_mammals.png'), p_mammals, width = 9, height = 11)
ggsave(here('figures/fig1_lms_mammals.pdf'), p_mammals, width = 9, height = 11)

# Fig 2 ----
plot_order = mod_coefs %>% 
  group_by(sp) %>% 
  slice(which.max(abs(estimate))) %>% 
  pull(bio_variable) %>% 
  table() %>% sort(decreasing = T) %>% names()

mod_coefs %>% 
  group_by(sp) %>% 
  slice(which.max(abs(estimate))) %>% 
  group_by(bio_variable) %>% 
  summarise(Positive = sum(estimate > 0, na.rm = T),
            Negative = sum(estimate < 0, na.rm = T)) %>% 
  gather("Direction", "count", Positive, Negative) %>% 
  mutate(bio_variable = factor(bio_variable, levels = plot_order)) %>% 
  ggplot(aes(x = bio_variable, y = count, fill = Direction)) +
  geom_bar(stat = 'identity', width = 0.4) +
  labs(y = 'Number of species', x = 'Variable that has the strongest effect') +
  coord_flip() +
  scale_fill_manual(values = c('orange', 'blue')) +
  theme(legend.position = c(0.7, 0.8))
ggsave(here('figures/fig2_strongest_est.png'), width = 6, height = 6)
ggsave(here('figures/fig2_strongest_est.pdf'), width = 6, height = 6)


# Fig partial R2
ggplot(par_r2, aes(x = var, y = partial_r2)) +
  ggparl::geom_boxjitter(outlier.color = NA, width = 0.5,
                         jitter.alpha = 0.6, jitter.shape = 21,
                         errorbar.draw = T, errorbar.length = 0.2,
                         fill = 'skyblue'
                         ) +
  labs(y = expression(paste("Partial adjusted ", R^{2})), 
       x = 'Climatic predictors') +
  theme(legend.position = 'none')
ggsave(here('figures/fig3_partial_r2.png'), width = 6, height = 5)
ggsave(here('figures/fig3_partial_r2.pdf'), width = 6, height = 5)

# how many species have the highest partial R2 of Temp? 60
group_by(par_r2, sp) %>% 
  slice(which.max(partial_r2)) %>% 
  pull(var) %>% 
  table() %>% 
  sort()

ggplot(par_r2, aes(y = partial_r2, x = sp)) +
  geom_line(aes(group = var, color = var)) +
  # geom_point(aes( color = var)) +
  theme(axis.text.x = element_blank()) +
  scale_y_log10()
