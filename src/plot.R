# Fig 1 ----
# mod_ceof2 from 03_lm.R
p_a = ggplot(mod_coefs, aes(x = bio_variable, y = sp)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  scale_fill_gradient2() +
  geom_text(aes(label = sig), nudge_y = -0.55) +
  labs(x = 'Climatic predictors', y = '', fill = 'Coefs') +
  theme(legend.position = c(0.44, 0.5),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size = 4),
        axis.ticks.y.left = element_line(size = 0.2),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))

p_b = ggplot(mod_coefs, aes(y = sp, x = adj.r.squared)) +
  geom_segment(aes(yend = sp), xend = 0, color = 'gray') +
  geom_point(color = 'blue', size = 0.8) +
  labs(x = expression(paste("Adjust ", R^{2})), 
       y = '') +
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
p = plot_grid(p_a, p_b, align = 'h', rel_widths = c(5, 1))
ggsave(here('figures/fig1_lms.pdf'), p, width = 8, height = 10)

# Fig 2 ----
mod_coefs %>% 
  group_by(sp) %>% 
  slice(which.max(abs(estimate))) %>% 
  group_by(bio_variable) %>% 
  summarise(sig_pos = sum(estimate > 0 & p.value < 0.05, na.rm = T),
            sig_neg = sum(estimate < 0 & p.value < 0.05, na.rm = T)) %>% 
  gather("sig", "count", sig_pos, sig_neg) %>% 
  mutate(bio_variable = factor(bio_variable, levels = c('Temp', 'Precip S.', 'Precip', 'Temp S.'))) %>% 
  ggplot(aes(x = bio_variable, y = count, fill = sig)) +
  geom_bar(stat = 'identity', width = 0.5) +
  labs(x = '', y = '') +
  coord_flip() +
  scale_fill_manual(values = c('orange', 'blue'))
ggsave(here('figures/fig2_strongest_est.pdf'), width = 6, height = 6)


# Fig partial R2
ggplot(par_r2, aes(x = var, y = partial_r2)) +
  ggparl::geom_boxjitter(outlier.color = NA, width = 0.5,
                         jitter.alpha = 0.6, jitter.shape = 21,
                         errorbar.draw = T, errorbar.length = 0.2,
                         fill = 'skyblue') +
  labs(y = expression(paste("Partial adjusted ", R^{2})), 
       x = 'Climatic predictors')
ggsave(heer('figures/fig3_partial_r2.pdf'), width = 6, height = 5)
