names(d)
summary(d$long) 
summary(d$lat)
# plot(d$mass_g, d$length_mm) # look okay now
hist(d$mass_g) # highly skewed
hist(d$mass_g_log10) # much better

# I'd round lat and long to 2 decimal (~1km).
# d = mutate(d, lat = round(lat, 2), long = round(long, 2))
# climate = group_by(d, lat, long) %>% 
#   summarise_at(.vars = vars(bio1, bio12, bio4, bio15), .funs = mean, na.rm = T) %>% 
#   ungroup()
# 
# pairs(climate)
# dat = select(d, lat, long, year, month, sp, mass_g, sex, seasonality) %>% 
#   left_join(climate, by = c('lat', 'long'))

p = ggplot(d, aes(x = bio1, y = mass_g_log10, color = year, group = year)) +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(~sp, ncol = 10, scales = 'free') +
  theme(legend.position = 'bottom') +
  labs(title = 'So idiosyncratic, not a good idea to look by year for each species; \n and most years only have 1 record (not shown).')
ggsave('figures/mass_temp_by_sp_yr.pdf', width = 20, height = 30)
# does not make sense to use yr as a predictor given the idiosyncratic slopes...


# pre-1980 and post-1980? seems more reasonable to use this
p2 = ggplot(d, aes(x = bio1, y = mass_g_log10, color = yr_80, group = yr_80)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(~sp, ncol = 10, scales = 'free') +
  theme(legend.position = 'bottom')
ggsave('figures/mass_temp_by_sp_1980.pdf', width = 20, height = 30)


# does the collection season matter for species?
season_aov = group_by(d, sp) %>% 
  do(broom::tidy(car::Anova(lm(mass_g_log10 ~ season, data = .))))
filter(season_aov, term == 'season') %>%
  mutate(sig = p.value < 0.05) %>% 
  pull(sig) %>% mean() 
# about 70% species' mass differ between seasons; so yes

p3 = ggplot(d, aes(x = season, y = mass_g)) +
  geom_jitter(width = 0.05, alpha = 0.2) +
  geom_boxplot() +
  facet_wrap(~sp, ncol = 10, scales = 'free') +
  theme(legend.position = 'bottom')
ggsave('figures/mass_by_season.pdf', p3, width = 20, height = 30)

p3 = ggplot(d, aes(x = yr_80, y = mass_g)) +
  geom_jitter(width = 0.05, alpha = 0.2) +
  geom_boxplot() +
  facet_wrap(~sp, ncol = 10, scales = 'free') +
  theme(legend.position = 'bottom')
ggsave('figures/mass_by_yr_1980.pdf', p3, width = 20, height = 30)

# how many records each file/sp has?
d %>% 
  group_by(sp, yr_80, lat, long) %>% 
  tally() %>% 
  arrange(desc(n))
# some species have few records pre-/post-1980

filter(d, sp == 'Molothrus_ater') %>% View()

x = d %>% 
  group_by(lat, long) %>% 
  summarise(nyr = n_distinct(year)) %>% 
  arrange(desc(nyr)) %>% 
  filter(nyr >= 10) %>% 
  ungroup() %>% 
  left_join(d) %>% 
  mutate(loc = paste(lat, long, sep = "_")) 
p = ggplot(x, aes(x = year, y = bio1)) +
  geom_point() + geom_smooth(method = 'lm') +
  facet_wrap(~loc, scales = 'free', ncol = 6) +
  labs(title = 'For locations with >= 10 year records, what are their temperature trends?')
ggsave('figures/temp_by_yr.pdf', p, width = 15, height = 30)
# hum some points, e.g. #3, #7, around LA, have decreased temperature over time??

d %>% select(long, lat) %>% unique() %>% nrow() # 26,739 unique lat/long combination
n_distinct(d$sp) # 175 species

d2 = d %>% 
  group_by(sp, long, lat) %>% 
  count() %>% 
  ungroup()
dim(d2) # 43893,  4

d2 %>% 
  group_by(sp) %>% 
  count() %>% arrange(nn)
# minimal collection location: 50
# how far away are they?

source('Rcode/map.R')
p4 = map_america +
  geom_point(data = d,
             aes(x = long, y = lat, color = yr_80), 
             inherit.aes = F, size = 0.2) +
  facet_wrap(~sp, ncol = 10) +
  theme(legend.position = 'bottom')
ggsave('figures/location_by_sp.pdf', p4, width = 15, height = 30)
# lost of species have limited geographic distribution, which justifies the use of LMM.


# z-standardize climatic variables
d[, c('bio1', 'bio12', 'bio4', 'bio15')] = scale(d[, c('bio1', 'bio12', 'bio4', 'bio15')])

cor(d[, c('bio1', 'bio12', 'bio4', 'bio15')]) # seems okay

# (1|sp/season): collection season nested within species to account for seasonal var in body size of each sp
lmm1 = lmer(mass_g_log10 ~ bio1 + (1|sp/season) + (bio1|sp), data = d, REML = F)
lmm2 = lmer(mass_g_log10 ~ bio1 + (1|sp/season) + (0 + bio1|sp), data = d, REML = F)
anova(lmm1, lmm2) # is correlation between random intercept and slope strong? No.
summary(lmm1)
summary(lmm2)
fixef(lmm2)

mass_temp_by_sp_lm_lmm_plot = ggplot(d, aes(x = bio1, y = mass_g_log10, color = sp)) +
  # geom_point(alpha = 0.3) +
  geom_smooth(method = 'lm', se = F, size = 0.5) +
  theme(legend.position = 'null') +
  labs(x = 'Annual average temperature',
       y = 'Log10 body mass in g') +
  geom_abline(intercept = fixef(lmm2)[1], slope = fixef(lmm2)[2], size = 1)
ggsave('figures/mass_temp_by_sp_lm_lmm_plot.pdf', plot = mass_temp_by_sp_lm_lmm_plot,
       width = 5, height = 5)

lmm3 = lmer(mass_g_log10 ~ bio1 * yr_80 + (1|sp/season) + (0 + bio1:yr_80|sp), data = d, REML = F)
summary(lmm3)
# no strong interaction between bio1 and yr_80, suggesting that mass ~ bio relationship did not change pre-/post-1980
# but body mass seems declined 

lmm4 = lmer(mass_g_log10 ~ yr_80 + bio1 + (1|sp/season) + (0 + bio1|sp), data = d, REML = F)
summary(lmm4)
# indeed, body mass decreased post-1980 significantly

# now add all predictors
system.time(lmm5 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio4 + bio12 + bio15 + 
                           (1 | sp/season) + (0 + yr_80 | sp) + (0 + bio1 | sp) + (0 + bio4 | sp) + 
                           (0 + bio12 | sp) + (0 + bio15 | sp), 
                         data = d, REML = F))
summary(lmm5)
MuMIn::r.squaredGLMM(lmm5)
rr2::R2(lmm5, resid = F, pred = F)

# remove bio15?
lmm6 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio4 + bio12 + # bio15 + 
               (1 | sp/season) + (0 + bio1 | sp) + (0 + bio4 | sp) + 
               (0 + bio12 | sp), 
             data = d, REML = F)
anova(lmm6, lmm5) # keep it

# does mass ~ bio1 relationship changed?
lmm7 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio1:yr_80 + bio4 + bio12 + bio15 + 
               (1 | sp/season) + (0 + bio1 | sp) + (0 + bio4 | sp) + 
               (0 + bio12 | sp) + (0 + bio15 | sp) + (0 + bio1:yr_80 | sp), 
             data = d, REML = F)
summary(lmm7)
# again, no strong interaction, suggesting no changes in mass ~ bio1 relationship

# correlation among random terms?
lmm8 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio4 + bio12 + bio15 + 
               (1 | sp/season) + (1 + bio1 + bio4 + bio12 + bio15 | sp), 
             data = d, REML = F)
summary(lmm8) # no strong correlations
anova(lmm5, lmm8) # no need

lmm5_no_bio1 <- lmer(mass_g_log10 ~ yr_80 + bio4 + bio12 + bio15 + 
                       (1 | sp/season) + (0 + bio4 | sp) + 
                       (0 + bio12 | sp) + (0 + bio15 | sp), 
                     data = d, REML = F)

lmm5_no_bio4 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio12 + bio15 + 
                       (1 | sp/season) + (0 + bio1 | sp) +
                       (0 + bio12 | sp) + (0 + bio15 | sp), 
                     data = d, REML = F)

lmm5_no_bio12 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio4 + bio15 + 
                        (1 | sp/season) + (0 + bio1 | sp) + (0 + bio4 | sp) + 
                        (0 + bio15 | sp), 
                      data = d, REML = F)

rr2::R2(mod = lmm5, mod.r = lmm5_no_bio1) # contribution of bio1
rr2::R2(mod = lmm5, mod.r = lmm5_no_bio4) # contribution of bio4
rr2::R2(mod = lmm5, mod.r = lmm5_no_bio12) # contribution of bio12

system.time(lmm9 <- lmer(mass_g_log10 ~ yr_80 + bio1 + bio4 + bio12 + bio15 + bio1:bio12 +  
                           (1 | sp/season) + (0 + bio1 | sp) + (0 + bio4 | sp) + 
                           (0 + bio12 | sp) + (0 + bio15 | sp), 
                         data = d, REML = F))
summary(lmm9)

niche = readxl::read_xlsx('data/NicheParameters175Species.xlsx')
niche = niche[1:175,]
names(niche) = c('sp', 'raw_records', 'used_records', 'mean_mass',
                 'mean_temp', 'breadth_temp', 'mean_precip', 'breadth_precip',
                 'mean_temp_var', 'breadth_temp_var', 'mean_precip_var', 'breadth_precip_var')
niche$sp = gsub(' ', '_', niche$sp)

plot(niche$mean_temp, niche$breadth_temp)

# because these niche parameters are all based on variables that already in the LMM,
# I am not sure about adding them into the LMM and look at their interactions with bio1.
coef_lmm5 = coef(lmm5)
niche = left_join(niche, rownames_to_column(coef_lmm5$sp, 'sp'))
select(niche, mean_temp:breadth_precip_var, bio1:bio15) %>% 
  pairs()
