# Packages ----

if(!require(xfun)) install.packages("xfun")
xfun::pkg_attach2(c("vegan", "tidyverse", "lme4", "rr2", "olsrr"))

# Functions ----

#' Generates cleaned body mass data files.  
#' 
#' Filters out records with spurious temperature and/or body size data.
#'
#' @param DataDir The directory with individual species data files.
#' @param species_path File that lists the species that we want to analyze.
#' @param output_dir Where to save the cleaned CSV files.
#' @param FromRecord Start of the row range to analyze; (integer).
#' @param ToRecord End of the row range to analyze; (integer, inclusive).
#'
generate_cleaned_CSV_files = function(DataDir, species_path, output_dir, FromRecord, ToRecord)
{
  #	t1 = read.csv(SpListFile, header=TRUE)
  
  RangeTbl = read.csv(species_path, header=TRUE, stringsAsFactors = FALSE)
  for (i in FromRecord:ToRecord)
  {
    SpName = trimws(as.character(RangeTbl[i,1]), which= "both")
    InpFileName = paste0("Env_", SpName, ".csv")
    inp_path = paste0(DataDir, InpFileName)
    
    print(paste(SpName, InpFileName, sep = ":"))
    
    Tbl2 = read.csv(inp_path, header=TRUE, stringsAsFactors=FALSE, fileEncoding="UTF-8" )
    print(dim(Tbl2))
    
    ## Remove data rows with temperature is too low. Specially to remove data
    ## with -9999 values in bio1
    NoDataIndices = which(Tbl2[['bio1']] <= -50)
    print(paste("Total -9999 rows ", length(NoDataIndices)))
    if (length(NoDataIndices) > 0 )
    {
      Tbl2 = Tbl2[-NoDataIndices,]
    }
    
    ## Find the truncated mean 10% of the data from upper and lower range removed. 
    ## Get the lower and upper bounds from the data, so that SD can be calculated. 
    Tbl2 = Tbl2[order(Tbl2$massing), ]
    trim_pct = 0.1
    TrimmedMean = mean(Tbl2$massing, na.rm = TRUE, trim = trim_pct)
    
    ## Changing the Midpoint to trimmed mean
    MidPoint = TrimmedMean		
    MinMass = MidPoint * 0.25
    MaxMass = MidPoint * 2.0
    
    ValidIndices = which(Tbl2[['massing']] >= MinMass & Tbl2[['massing']] <= MaxMass)
    Tbl3 = Tbl2[ValidIndices,]
    Tbl3 = Tbl3[Tbl3$seasonality != '',]
    
    print(dim(Tbl3))
    
    write.csv(Tbl3, paste0(output_dir, SpName, '.csv'), row.names=FALSE)
  }
}

#' Linear regression after z-sclae predictors
#' @param x species name
#' @param fm formula character string
#' @param df a data frame with the following columns: bio1, bio12, bio4, bio15, season
#' @return a data frame with the final model, coefs, and partial R2s. etc.
lm_1_sp = function(x = "Accipiter_cooperii", 
                   fm = "mass_g_log10 ~  bio1 + bio12 + bio4 + bio15 + season", 
                   df = d)
{
  cat(x, "\n")
  df = filter(df, sp == x)
  
  # vif and conditional indices; conditional indices requires un-centered data...
  mod_0 = lm(as.formula(fm), data = df)
  mod_diag = ols_coll_diag(mod_0)
  coll_prob = any(mod_diag[[2]]$`Condition Index` > 30)
  # which variables cause multicollinearity?
  if(coll_prob){
    ci = filter(mod_diag$eig_cindex, `Condition Index` > 30) %>% 
      select(intercept, starts_with("bio")) %>% 
      as.matrix()
    which_coll = colnames(ci)[ci[1,] > 0.5]
  } else {
    which_coll = NA
  }
  
  # output data frame
  out = tibble(
    n_row = nrow(df),
    vif_raw = mod_diag[1],
    cond_raw = mod_diag[2],
    multicol_prob = coll_prob,
    which_coll = list(which_coll)
  )
  
  # scale predictors
  df0 = df
  df[, c("bio1", "bio12", "bio4", "bio15")] = scale(df[, c("bio1", "bio12", "bio4", "bio15")]) 
  # fit model again
  if(coll_prob && # has multicol problem
     (!identical(which_coll, c("intercept", "bio4")) & # except intercept and bio4 coll
     any(grepl("^bio", which_coll)) # except coll between intercept (Spring) and season; i.e. must has bioX variable
     )){ 
    # remove bio4
    mod = lm(mass_g_log10 ~  bio1 + bio12 + bio15 + season, data = df)
    mod_1 = lm(mass_g_log10 ~  bio1 + bio12 + bio15 + season, data = df0) # raw data
    coll_d2 = ols_coll_diag(mod_1)
    coll_remain = any(coll_d2$eig_cindex$`Condition Index` > 30)
  } else {
    mod = lm(as.formula(fm), data = df)
    coll_remain = FALSE
  }
  
  mod_coef <- broom::tidy(mod)
  
  # output data frame
  out = mutate(out, coll_remain = coll_remain,
         final_mod = list(mod), 
         coef = list(mod_coef), # coef and p values of predictors
         partial_r2_bio1 = NA,
         partial_r2_bio4 = NA,
         partial_r2_bio12 = NA,
         partial_r2_bio15 = NA
         )
  
  # partial R2 of predictors
  bio_vars = grep("^bio", x = all.vars(formula(mod)), value = T)
  for(i in bio_vars){
    new_fm = paste0(". ~ . - ", i)
    mod.r = update(mod, new_fm)
    out[, paste0("partial_r2_", i)] = rr2::partialR2adj(mod = mod, mod.r = mod.r)$R2.adj
  }
  out = bind_cols(out, broom::glance(mod)) %>% 
    bind_cols(filter(mod_coef, term %in% paste0("bio", c(1, 4, 12, 15))) %>% 
    select(-std.error, -statistic) %>% 
    gather('var', 'value', estimate, p.value) %>% 
    unite('coefs', term, var) %>% 
    spread('coefs', 'value'))
  
  out
}
