# Packages ----

if(!require(xfun)) install.packages("xfun")
if(!require(rr2)) devtools::install_github("arives/rr2")
xfun::pkg_attach2(c("vegan", "tidyverse", "lme4", "rr2"))

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
  cat("sp = ", x, "\t")
  df = filter(df, sp == x)
  # scale predictors
  df[, c("bio1", "bio12", "bio4", "bio15")] = scale(df[, c("bio1", "bio12", "bio4", "bio15")]) 
  # df = rename(df, temp = bio1, precip = bio12, 
  #             temp_seasonality = bio4, precip_seasonality = bio15)
  
  mod = lm(as.formula(fm), data = df)
  # VIF for multicollinearity
  gvif = car::vif(mod)[, "GVIF"]
  reduced = any(gvif > 3)
  while(any(gvif > 3)){
    # remove the predictor has the highest VIF from the formula
    fm = gsub(paste0(names(which.max(gvif)), " [+] "), "", fm)
    # refit the model
    mod = lm(as.formula(fm), data = df)
    # VIF again
    gvif = car::vif(mod2)[, "GVIF"]
  }
  
  # output data frame
  out = tibble(
               n_row = nrow(df),
               multicol_prob = reduced,
               full_mod = list(mod), 
               coef = list(broom::tidy(mod)), # coef and p values of predictors
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
  out = bind_cols(out, broom::glance(mod))
  
  out
}
