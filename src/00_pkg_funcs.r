#------------
# Copyright (C) 2018 Brian J. Stucky.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#------------


#------------
# Provides functions that implement various conditioning and collinearity
# diagnostics for linear models.  The methods are based on those developed in
# Belsley 1991.
#------------


#------
# Tests whether a matrix is singular.
#------
is_singular = function(matx) {
  return( class(try(solve(matx), silent=T)) == 'try-error' )
}


#------
# Scales a matrix so all columns have unit length under the 2-norm.
#------
scale_matrix = function(matx) {
  for (colnum in 1:ncol(matx)) {
    colvec = matx[,colnum]
    # If the norm of the vector is 0 (that is, the vector is a 0 vector), leave
    # it as is.
    if (sum(colvec) != 0) {
      matx[,colnum] = colvec / norm(colvec, type='2')
    }
  }
  
  return(matx)
}


#------
# Gets the condition number, condition indices, and variance decomposition
# proportions associated with a linear model formula and predictor variables.
#------
collin_diags = function(model_spec, data=environment(model_spec)) {
  d_mtx = scale_matrix(model.matrix(model_spec, data=data))
  
  xtx = t(d_mtx) %*% d_mtx
  if (is_singular(xtx)) {
    result = list(c_num=Inf, c_indices=NULL, VDPs=NULL)
    class(result) = c('CollinDiagsResults', class(result))
    return(result)
  }
  
  # Because we tested for the invertibility of X'X, X'X will not have 0 as an
  # eigenvalue, so all singular values of X will be nonzero.
  svdv = svd(d_mtx, nv=ncol(d_mtx))
  svals = svdv$d
  
  c_num = max(svals) / min(svals)
  
  # Calculate the condition indices.
  c_indices = max(svals) / svals
  
  # Calculate the variance decomposition proportions for each coefficient
  # in the regression model.
  var_decomp = matrix(nrow=0, ncol=ncol(d_mtx))
  for (varnum in 1:ncol(d_mtx)) {
    phis = (svdv$v[varnum,] ^ 2) / (svals ^ 2)
    var_decomp = rbind(var_decomp, phis / sum(phis))
  }
  colnames(var_decomp) = round(c_indices, 2)
  rownames(var_decomp) = colnames(d_mtx)
  
  result = list(
    c_num=c_num,
    c_indices=sort(c_indices, decreasing=T),
    VDPs=var_decomp
  )
  class(result) = c('CollinDiagsResults', class(result))
  
  return(result)
}

#------
# Defines a custom print function for collinearity diagnostics result objects.
#------
print.CollinDiagsResults = function(result) {
  cat('Condition number of design matrix: ', result$c_num, '\n\n', sep='')
  cat('Condition indices:\n')
  cat(result$c_indices, sep=', ')
  cat('\n\nVariance Decomposition Proportions:\n')
  print(round(result$VDPs, 2))
}


#------
# Generates cleaned body mass data files.  Filters out records with spurious
# temperature and/or body size data.
#
# DataDir (str): The directory with individual species data files.
# sp_table (dataframe): A data frame of metadata for the species that we want to analyze.
# output_dir: Where to save the cleaned CSV files.
#------
generateCleanedCSVFiles = function(DataDir, sp_table, output_dir) {
  for (i in 1:nrow(sp_table)) {
    SpName = trimws(as.character(sp_table[i,1]), which= "both")
    InpFileName = paste0("Env_", SpName, ".csv")
    inp_path = paste0(DataDir, InpFileName)
    
    print(paste0('Cleaning data for ', SpName, ' (', InpFileName, ')...'))
    
    Tbl2 = read.csv(inp_path, header=TRUE, stringsAsFactors=FALSE, fileEncoding="UTF-8" )
    startcnt = nrow(Tbl2)
    
    # Remove data rows with temperature is too low. Specially to
    # remove data with -9999 values in bio1
    NoDataIndices = which(Tbl2[['bio1']] <= -50)
    if (length(NoDataIndices) > 0 ) {
      Tbl2 = Tbl2[-NoDataIndices,]
    }
    
    
    # Find the truncated mean 10% of the data from upper and lower
    # range removed.  Get the lower and upper bounds from the data,
    # so that SD can be calculated. 
    Tbl2 = Tbl2[order(Tbl2$massing), ]
    trim_pct = 0.1
    TrimmedMean = mean(Tbl2$massing, na.rm = TRUE, trim = trim_pct)
    
    # Calculate the mass cutoff values.
    MidPoint = TrimmedMean        
    MinMass = MidPoint * 0.25
    MaxMass = MidPoint * 2.0
    
    ValidIndices = which(Tbl2[['massing']] >= MinMass & Tbl2[['massing']] <= MaxMass)
    Tbl3 = Tbl2[ValidIndices,]
    Tbl3 = Tbl3[Tbl3$seasonality != '',]
    
    print(paste0(startcnt - nrow(Tbl3), ' record(s) removed.'))
    
    write.csv(Tbl3, paste0(output_dir, SpName, '.csv'), row.names=FALSE)
  }
}


#------
# Writes information about collinearity in a model to a file.
#------
write_collin_report = function(fname, res) {
  if (res$c_num == Inf) {
    return()
  }
  
  sink(fname)
  print(res)
  sink()
}


#----
# Standardizes one or more columns in a dataframe.  Predictors are centered to
# have mean == 0 and scaled to have s.d. == 1.0. 
#----
standardize_cols = function(df, colnames) {
  for (colname in colnames) {
    df[[colname]] = as.vector(scale(df[[colname]], scale=TRUE, center=TRUE))
  }
  
  return(df)
}


#----
# Generates model fit diagnostic outputs.
#----
PlotResiduals = function(CurModel, SpName, mod_name, resDir) {
  jpegfileName = paste(resDir, SpName, "_", mod_name, "_res-qq.jpg", sep = "")
  jpeg(filename=jpegfileName, width=8, height=5, units = "in", res = 150)
  par(mfrow=c(1,2))
  resids = rstandard(CurModel)
  plot(CurModel$fitted.values, resids, ylab="Standardized residuals", xlab="Fitted values",
       main = paste(SpName, mod_name, sep ="_"))
  abline(h=0, col="red", lty=2, lwd=1.4)
  qqnorm(resids, main = paste(SpName, mod_name, sep ="_"))
  abline(a=0, b=1, col="red", lty=2, lwd=1.4)
  dev.off()
  
  colnames = names(CurModel$model)
  colnames = colnames[2:length(colnames)]
  
  jpegfileName = paste(resDir, SpName, "_", mod_name, "_res-preds.jpg", sep = "")
  jpeg(filename=jpegfileName, width=8, height=6, units = "in", res = 150)
  par(mfrow=c(ceiling(length(colnames) / 3),3))
  
  for (colname in colnames) {
    t1 = data.frame(predictor=CurModel$model[[colname]], residual=resids)
    
    t2 = complete.cases(t1)
    t1 = t1[t2,]
    plot(t1$predictor, t1$residual, xlab=colname, ylab="Standardized residuals",
         main = paste(SpName, mod_name, colname, sep ="_"))
    
    abline(h=0, col="red", lty=2, lwd=1.4)
  }
  
  dev.off()
  
  ## Data table with residual and fitted values. 
  mat1 = cbind(CurModel$res, CurModel$fitted.values, CurModel$model)
  resFileName = paste(resDir, SpName, "_", mod_name, "_res.csv", sep = "")
  print(resFileName)
  names(mat1)[1:2] = c("residuals", "fitted.values")
  write.csv(mat1, resFileName, row.names=FALSE)
}


#----
# Generates a matrix containing a summary of a fitted model.
#----
GetModelSummary = function(model, resDir, SpName, modelName) {
  ## Get the coefficients and store it in matrix
  OpMat1 = summary(model)$coefficients
  
  ## Get the string of significant variables in variable SigVar
  ## To add (+) or (-) sign in the significant variable. 
  pos = which(OpMat1[-1,4] < 0.05 & OpMat1[-1,1] > 0)
  neg = which(OpMat1[-1,4] < 0.05 & OpMat1[-1,1] < 0)
  SigVarNamePos = ''
  SigVarNameNeg = ''
  
  if (length(pos) > 0) {
    SigVarNamePos = paste(names(pos), '(+)', collapse = "; ", sep = "")
  }
  if (length(neg) > 0) {
    SigVarNameNeg = paste(names(neg), '(-)', collapse = "; ", sep = "")
  }
  
  SigVarName = paste(SigVarNamePos, SigVarNameNeg, sep = "; ")
  
  ## Add the row names as one of the column in the matrix
  OpMat1 = cbind(row.names((summary(model)$coefficients)), OpMat1)
  
  ## Add p values and r squares to the matrix
  pValue = pf(summary(model)$fstatistic[1L], summary(model)$fstatistic[2L], summary(model)$fstatistic[3L], lower.tail = FALSE)
  Rsquared = summary(model)$r.squared
  AdjRsquared = summary(model)$adj.r.squared
  deg_fr = model$df
  
  pval = c("pval", pValue, "", "", "")
  r_square = c("r-square",Rsquared, "", "", "")
  adj_r_square = c("adj.r.square", AdjRsquared, "", "", "")
  degree_freedom = c("degrees_of_freedom", deg_fr, "", "", "")
  
  OpMat1 = rbind(OpMat1, pval, r_square, adj_r_square, degree_freedom)
  row.names(OpMat1) = NULL
  
  ## save the models in results directory
  OpFileName = paste(resDir, SpName, "_", modelName, ".csv", sep = "")
  write.csv(OpMat1, OpFileName, row.names=FALSE)
  
  # Now get the last four rows from each model and transpose it to
  # columns and attach it to output matrix.
  TotRows= dim(OpMat1)[1]
  nt1 = OpMat1[(TotRows - 3):TotRows, 2]
  nt1 = c(nt1, SigVarName)
  
  return(nt1)
}


#------
# Runs model fits for all rows in a master species list table.
#
# DataDir (str): The directory with individual species data files.
# sp_table (dataframe): A data frame of metadata for the species that we want to analyze.
# output_dir (str): Directory for the modeling results.
# model_specs: A list of model formulas to fit.
# res_models: Model names for which to generate fit diagnostics.
#------
Model_All = function(DataDir, sp_table, output_dir, model_specs, res_models) {
  collin_dir = paste0(output_dir, 'collin_results/')
  if (!dir.exists(collin_dir)) { dir.create(collin_dir) }
  diags_dir = paste0(output_dir, 'diagnostics/')
  if (!dir.exists(diags_dir)) { dir.create(diags_dir) }
  modfits_dir = paste0(output_dir, 'model_fits/')
  if (!dir.exists(modfits_dir)) { dir.create(modfits_dir) }
  
  TotSps = dim(sp_table)[1]
  num_rescols = 1 + length(model_specs) * 7
  ResMat = matrix(0, nrow=0, ncol = num_rescols)
  CollinTbl = matrix("", nrow=0, ncol=13)
  preds_list = c(
    '(Intercept)', 'bio1', 'bio12', 'bio4', 'bio15', 'seasonalitySpring',
    'seasonalitySummer', 'seasonalityWinter'
  )
  
  for (i in 1:nrow(sp_table)) {
    SpName = trimws(as.character(sp_table[i,1]), which= "both")
    InpFileName = paste(DataDir, SpName, ".csv", sep = "")
    
    print(paste0('Analyzing data for ', SpName, ' (', InpFileName, ')...'))
    
    Tbl3 = read.csv(InpFileName, header=TRUE, stringsAsFactors=FALSE, fileEncoding="UTF-8" )
    
    # Initialize the output results row.
    CurRow = SpName
    
    # Do each model fit and get the fit statistics. 
    for (mod_name in names(model_specs)) {
      model_spec = model_specs[[mod_name]]
      
      collin_res = collin_diags(model_spec, Tbl3)
      c_num = collin_res$c_num
      
      ncols = ncol(collin_res$VDPs)
      newrow = rep(0, length(preds_list))
      row_exist = which(preds_list %in% row.names(collin_res$VDPs))
      newrow[row_exist] = collin_res$VDPs[,ncols]
      CollinTbl = rbind(
        CollinTbl,
        c(SpName, mod_name, as.character(collin_res$c_indices[1:3]), as.character(newrow))
      )
      
      if (c_num > 30) {
        fname = paste0(collin_dir, SpName, '-', mod_name, '.txt')
        write_collin_report(fname, collin_res)
      }
      
      # Standardardize the continuous variables before model fitting.
      Tbl3_stdzd = standardize_cols(Tbl3, c('bio1', 'bio12', 'bio4', 'bio15'))
      Tbl3_stdzd$massing = scale(Tbl3_stdzd$massing)
      mod = lm(model_spec, data=Tbl3_stdzd)
      
      # This is to plot the residuals and save the residuals and
      # fitted values in a text file. 
      if (length(which(mod_name %in% res_models)) > 0 ) {
        PlotResiduals(mod, SpName, mod_name, diags_dir)
      }
      mod_summary = GetModelSummary(mod, modfits_dir, SpName, mod_name)
      
      CurRow = c(CurRow, mod_summary, c_num, AIC(mod))
    }
    
    ResMat = rbind(ResMat, CurRow)
  }
  
  res_df = as.data.frame(ResMat)
  
  # Generate the column names for the results data frame.
  mod_colname_suffixes = c(
    "_pval", "_r_square", "_adj_r_square", "_degree_freedom", "_sig_vars", '_condnum', "_AIC"
  )
  colnames = c('Species')
  for (mod_name in names(model_specs)) {
    colnames = c(colnames, paste(mod_name, mod_colname_suffixes, sep=''))
  }
  names(res_df) = colnames
  write.csv(res_df, paste0(output_dir, 'model_summaries.csv'), row.names=FALSE)
  CollinTbl = as.data.frame(CollinTbl)
  names(CollinTbl) = c('SpeciesName', 'mod', 'ci_1', 'ci_2', 'ci_3', preds_list)
  write.csv(CollinTbl, paste(output_dir, 'condition_summaries.csv', sep=''), row.names=FALSE)
  
  return(ResMat)
}


# added by D Li ----

if(!require(xfun)) install.packages("xfun")
if(!require(ggparl)) devtools::install_github('erocoar/ggparl')
xfun::pkg_attach2(c("tidyverse", "rr2", "cowplot", "broom", "here", "vegan"))

#' Linear regression after z-sclae predictors
#' @param x species name
#' @param mod_type which model to use?
#' @param data_path where are the data?
#' @param mass_transf how to transfer mass?
#' @param verbose do you want to see which sp is running?
#' @return a data frame with the final model, coefs, and partial R2s. etc.
lm_1_sp = function(x = "Accipiter_cooperii", 
                   mod_type = 'm8', 
                   data_path = "data_output/cleaned_bodymass_data/",
                   mass_transf = c('log10', 'scale_', 'scale_log10'),
                   verbose = FALSE)
{
  if(verbose) cat(x, "\n")
  
  mass_transf = match.arg(mass_transf)
  df = read.csv(paste0(data_path, x, '.csv'), stringsAsFactors = F) %>% as.tibble()
  # rename ...
  df = rename(df, season = seasonality, mass_in_g = massing)
  df = mutate(df, season = factor(season, levels = c('Spring', 'Summer', 'Fall', 'Winter')))
  # scale predictors to have mean 0 and sd 1
  df[, c("bio1", "bio12", "bio4", "bio15")] = scale(df[, c("bio1", "bio12", "bio4", "bio15")]) 
  if(mass_transf == 'log10') df$Y = log10(df$mass_in_g)
  if(mass_transf == 'scale_') df$Y = scale(df$mass_in_g, center = FALSE)
  if(mass_transf == 'scale_log10') df$Y = scale(log10(df$mass_in_g), center = FALSE)
  # log_10 body mass
  
  if(mod_type == 'm8') fm = 'Y ~ bio1 + bio12 + bio4 + bio15 + season'
  if(mod_type == 'm10') fm = 'Y ~ bio1 + bio12 + bio15 + season'
  
  # fit the model
  mod = lm(as.formula(fm), data = df)
  mod_coef = broom::tidy(mod)
  mod_r2 = broom::glance(mod)
  
  # output data frame
  out = tibble(
    n_row = nrow(df),
    final_mod = list(mod), 
    coef = list(mod_coef),
    partial_r2_bio1 = NA,
    partial_r2_bio4 = NA,
    partial_r2_bio12 = NA,
    partial_r2_bio15 = NA,
    partial_r2_bio1_bio4 = NA,
    partial_r2_bio12_bio15 = NA
  ) %>% 
    bind_cols(mod_r2)
  
  # partial R2 of predictors
  bio_vars = grep("^bio", x = all.vars(formula(mod)), value = T)
  for(i in bio_vars){
    new_fm = paste0(". ~ . - ", i)
    mod.r = update(mod, new_fm)
    out[, paste0("partial_r2_", i)] = rr2::partialR2adj(mod = mod, mod.r = mod.r)$R2.adj
  }
  if(all(c('bio1', 'bio4') %in% bio_vars)){
    mod.r = update(mod, ". ~ . - bio1 - bio4")
    out[, "partial_r2_bio1_bio4"] = rr2::partialR2adj(mod = mod, mod.r = mod.r)$R2.adj
  }
  if(all(c('bio12', 'bio15') %in% bio_vars)){
    mod.r = update(mod, ". ~ . - bio12 - bio15")
    out[, "partial_r2_bio12_bio15"] = rr2::partialR2adj(mod = mod, mod.r = mod.r)$R2.adj
  }
  
  # extract coef and p-values
  mod_coef2 = filter(mod_coef, term %in% paste0("bio", c(1, 4, 12, 15))) %>% 
    select(-std.error, -statistic) %>% 
    gather('var', 'value', estimate, p.value) %>% 
    unite('coefs', term, var) %>% 
    spread('coefs', 'value')
  
  out = bind_cols(out, mod_coef2)
  
  out
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, color.threshold = 0.5, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  z = na.omit(data.frame(x = x, y = y))
  r <- cor(z$x, z$y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) 
    cex.cor <- 0.9/strwidth(txt)
  color.cor = ifelse(abs(r) > color.threshold, "red", "black")
  text(0.5, 0.5, txt, cex = 2, col = color.cor)
  # text(0.5, 0.5, txt, cex = 7*cex.cor * r)
}

# Panel of pairs() function
#' \code{panel.hist} to plot absolute value of correlations
#' 
#' @export
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
}
