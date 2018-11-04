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
sum(lm_by_sp$multicol_prob) # 97 sp has potential multicollinearity problem
sum(lm_by_sp$coll_remain) # 0 sp has potential multicollinearity problem after removing bio4

hist(lm_by_sp$adj.r.squared)

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
