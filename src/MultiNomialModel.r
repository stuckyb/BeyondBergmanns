## Pass the table with variables, 
## raw_dt = dt1[, c("Species", "es_bio1", "pval_bio1", "adjPval_bio1", "sig_bio1", "Mean_used", "bio1_mean", "bio1_breadth")]

library(nnet)
library(lmtest)


## basically, multinom function does not take data as passed in variable. I do not understand why. So 

MetaFileName = "/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv"  
confintval = 0.95  
TypeOfPVal = "FDR"


MultiMod <- function(MetaFileName, confintval, TypeOfPVal, VarNameList, years)
{
  #dt1 = read.csv("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", header = TRUE)
  dt1 = read.csv(MetaFileName, header = TRUE) 
  print(names(dt1))
  if (TypeOfPVal == "FDR")
  {
    #VarList = list(
    #  'b1' = c("SpName", "Mean_used", "es_bio1", "pval_bio1", "adjPval_bio1", "sig_bio1_f", "sig_bio1_nf", "bio1_mean", "bio1_breadth"),
    #  'b2' = c("SpName", "Mean_used", "es_bio12", "pval_bio12", "adjPval_bio12", "sig_bio12_f", "sig_bio12_nf", "bio12_mean", "bio12_breadth"),
    #  'b3' = c("SpName", "Mean_used", "es_bio4", "pval_bio4", "adjPval_bio4", "sig_bio4_f", "sig_bio4_nf", "bio4_mean", "bio4_breadth"),
    #  'b4' = c("SpName", "Mean_used", "es_bio15", "pval_bio15", "adjPval_bio15", "sig_bio15_f", "sig_bio15_nf", "bio15_mean", "bio15_breadth")
    #)
    VarList = list(
      'b1' = c("SpName", "Mean_used", paste("es_bio1", years, sep = ""), paste("pval_bio1", years, sep = ""), paste("adjPval_bio1", years, sep = ""), paste("sig_bio1", years, "_f", sep = ""), 
	         paste("sig_bio1", years, "_nf", sep =""), paste("bio1", years, "_mean", sep =""), paste("bio1", years, "_breadth", sep = "")),
      'b2' = c("SpName", "Mean_used", paste("es_bio12", years, sep=""), paste("pval_bio12", years, sep=""), paste("adjPval_bio12", years, sep=""), paste("sig_bio12", years, "_f", sep =""), 
	         paste("sig_bio12", years, "_nf", sep = ""), paste("bio12", years, "_mean", sep =""), paste("bio12", years, "_breadth", sep = "")),
      'b3' = c("SpName", "Mean_used", paste("es_bio4", years,  sep = "" ), paste("pval_bio4", years, sep = ""), paste("adjPval_bio4", years, sep = ""), paste("sig_bio4", years, "_f", sep = ""), 
	         paste("sig_bio4", years, "_nf", sep = ""), paste("bio4", years, "_mean", sep = ""), paste("bio4", years, "_breadth", sep = "")),
      'b4' = c("SpName", "Mean_used", paste("es_bio15", years, sep = ""), paste("pval_bio15", years, sep = ""), paste("adjPval_bio15", years, sep = ""), paste("sig_bio15", years, "_f", sep = ""), 
	         paste("sig_bio15", years, "_nf", sep = ""), paste("bio15", years, "_mean", sep = ""), paste("bio15", years, "_breadth", sep = ""))
    )

  } else
  {
    VarList = list(
      'b1' = c("SpName", "Mean_used", "es_bio1", "pval_bio1", "adjPval_bio1", "sig_bio1_nf", "sig_bio1_f", "bio1_mean", "bio1_breadth"),
      'b2' = c("SpName", "Mean_used", "es_bio12", "pval_bio12", "adjPval_bio12", "sig_bio12_nf", "sig_bio12_f", "bio12_mean", "bio12_breadth"),
      'b3' = c("SpName", "Mean_used", "es_bio4", "pval_bio4", "adjPval_bio4", "sig_bio4_nf", "sig_bio4_f", "bio4_mean", "bio4_breadth"),
      'b4' = c("SpName", "Mean_used", "es_bio15", "pval_bio15", "adjPval_bio15", "sig_bio15_nf", "sig_bio15_f", "bio15_mean", "bio15_breadth")
    )
  }
  
  
  #ModelList = list(
  #  'm1' = levelDtNS ~ Mean_used + bio1_mean + bio1_breadth,
  #  'm2' = levelDtNS ~ Mean_used + bio12_mean + bio12_breadth,
  #  'm3' = levelDtNS ~ Mean_used + bio4_mean + bio4_breadth,
  #  'm4' = levelDtNS ~ Mean_used + bio15_mean + bio15_breadth
  #)
  
  
  print("Going in switch")
  
  switch(years,

  "_2yr" = {
	print("Going in 1yr")
    ModelList = list(
    'm1' = levelDtNS ~ Mean_used + bio1_2yr_mean + bio1_2yr_breadth,
    'm2' = levelDtNS ~ Mean_used + bio12_2yr_mean + bio12_2yr_breadth,
    'm3' = levelDtNS ~ Mean_used + bio4_2yr_mean + bio4_2yr_breadth,
    'm4' = levelDtNS ~ Mean_used + bio15_2yr_mean + bio15_2yr_breadth)
    },
  "_5yr" = {
    ModelList = list(
    'm1' = levelDtNS ~ Mean_used + bio1_5yr_mean + bio1_5yr_breadth,
    'm2' = levelDtNS ~ Mean_used + bio12_5yr_mean + bio12_5yr_breadth,
    'm3' = levelDtNS ~ Mean_used + bio4_5yr_mean + bio4_5yr_breadth,
    'm4' = levelDtNS ~ Mean_used + bio15_5yr_mean + bio15_5yr_breadth)
	},
	{
	ModelList = list(
    'm1' = levelDtNS ~ Mean_used + bio1_mean + bio1_breadth,
    'm2' = levelDtNS ~ Mean_used + bio12_mean + bio12_breadth,
    'm3' = levelDtNS ~ Mean_used + bio4_mean + bio4_breadth,
    'm4' = levelDtNS ~ Mean_used + bio15_mean + bio15_breadth)
	}
	
  )  ## switch
  
	print(ModelList)


  #VarNameList = c("bio1", "bio12", "bio4", "bio15")
  # raw_dt = dt1[, c("Species", "Mean_used", "es_bio1", "pval_bio1", "adjPval_bio1", "sig_bio1", "sig_bio1_nf", "bio1_mean", "bio1_breadth")]
  
  Op = matrix(0, nrow=0, ncol = 23)
  ## op_mc = Output data frame for model comparison
  Op_mc = matrix("", nrow=length(VarNameList), ncol = 6)
  ## op_pc = Output data frame for predictor comparison
  Op_pc = matrix("", nrow=0, ncol = 7)
  for (count in 1:length(VarList))
  {
    print(count)
    print(VarList[[count]])
    
    raw_dt = dt1[, VarList[[count]]]
    print(names(raw_dt))
    
    ## As per the new method remove the predictors which are not contributing in the model. This will be denoted by "Exl" meaning excluded. (March 2018)
    ## remove all the rows with sig_<varName>_<f/nf> with value "Exl"
    ## whichCol
    exl = which(raw_dt[,6] == "Exl")
    ## Sixth column will categorised classificance of p values. and the value "Exl" suggests that this row is excluded for further analysis. 
    if (length(exl) > 0)
    {
      raw_dt = raw_dt[-exl,]
    }
    print(dim(raw_dt))
    
    
    
    ## For non significant
    #raw_dt$levelDtNS = relevel(droplevels(raw_dt[,6]), ref = "NS")

    ## For -ve significance
    raw_dt$levelDtNS = relevel(droplevels(raw_dt[,6]), ref = "Neg")
    #model1 = multinom(levelDtNS ~ Mean_used + bio1_mean + bio1_breadth, raw_dt)
    #mcount = (count - 1) * 3
    model1 = multinom(ModelList[[count]], raw_dt)
    mod_null = multinom(levelDtNS ~ 1, raw_dt)
    mod_nullAIC = summary(mod_null)$AIC
    pseudo_rsquare = (logLik(mod_null) - logLik(model1)) / logLik(mod_null)
    lrt_res = lrtest(model1, mod_null)
    chisqval = lrt_res$Chisq[2]
    chipval = lrt_res$'Pr(>Chisq)'[2]
    
    ## exploring the contribution of individual predictors in the model
    fterms = attr(formula(model1), 'term.labels')
    pc_pvals = c()
    pc_chisqs = c()
    for (fterm in fterms) {
      drop1_model = update(model1, as.formula(paste0('~ . - ', fterm)))
      lrt_res = lrtest(model1, drop1_model)
      pc_chisqs = c(pc_chisqs, lrt_res$Chisq[2])
      pc_pvals = c(pc_pvals, lrt_res$'Pr(>Chisq)'[2])
      print(drop1_model)
      print(lrt_res)
    }
    Op_pc = rbind(Op_pc, c(VarNameList[[count]], pc_pvals, pc_chisqs))
    
    ## For Non significant
    ## Op1 = MultiLogistic1(model1, confintval, "NS", VarNameList[count])
    ## for reference as negative significance
    Op1 = MultiLogistic1(model1, confintval, "Neg", VarNameList[count])
    Op_mc[count,] = c(VarNameList[count], summary(model1)$AIC, mod_nullAIC, pseudo_rsquare, chisqval, chipval)
    Op = rbind(Op, Op1)		
    ## readline()
    
    ## plot multinomial model plot
    OpFileName = paste("./species_data/Multinom_nfdr", VarNameList[count], ".png", sep = "")
    #OpFileName = paste("./species_data/Multinom_", VarNameList[count], ".pdf", sep = "")
    prednames = names(raw_dt)[c(2,8,9)]
    # Multinom_plot(prednames, model1, OpFileName, raw_dt)
    
  }
  ## To get the names of variables, 
  ## varNames = dimnames(exp(coef(model1)))[[2]]
  varNames = c("Intercept", "AvgBodyMass", "VarMean", "VarBreadth")
  ColNames  = c("Variable", "Referennce", "ComparedTo", paste("coef", varNames, sep = "_"), paste("LB", varNames, sep = "_"), paste("UB", varNames, sep = "_"), paste("Pval", varNames, sep = "_"), "ConfidenceVal", "Deviance", "ModelAIC", "DegFreedom")
  Op_df = as.data.frame(Op)
  names(Op_df)  = ColNames
  
  Op_mc_df = as.data.frame(Op_mc)
  names(Op_mc_df) = c("VarName", "AIC_full", "AIC_null", "pseudo_r2", "Chisq", "pval")
  
  Op_pc_df = as.data.frame(Op_pc)
  names(Op_pc_df) = c("VarName", "pval_bodymass", "pval_centroid", "pval_breadth", "chiq_bodymass", "chiq_centroid", "chiq_breadth")
  return(list(Op_df, Op_mc_df, Op_pc_df))
  
}


## MultiMod (MetaFileName, confintval, TypeOfPVal, VarNameList, years)


## For 1 year
## VarNameList = c("bio1", "bio12", "bio4", "bio15")
## tt2 = MultiMod("meta_model_data_1yr_sex.csv", 0.95, "FDR", VarNameList, "")
## write.csv(tt2[[1]], "Multinomial_Summary-Neg_ref_1yr.csv", row.names=FALSE)
## write.csv(tt2[[2]], "Multinomial_MC-Neg_ref_1yr.csv", row.names=FALSE)

## For 2 year
## VarNameList = c("bio1_2yr", "bio12_2yr", "bio4_2yr", "bio15_2yr")
## tt2 = MultiMod("meta_model_data_2yr_sex.csv", 0.95, "FDR", VarNameList, "_2yr")
## write.csv(tt2[[1]], "Multinomial_Summary-Neg_ref_2yr.csv", row.names=FALSE)
## write.csv(tt2[[2]], "Multinomial_MC-Neg_ref_2yr.csv", row.names=FALSE)


## For 5 year
## VarNameList = c("bio1_5yr", "bio12_5yr", "bio4_5yr", "bio15_5yr")
## tt2 = MultiMod("meta_model_data_5yr_sex.csv", 0.95, "FDR", VarNameList, "_5yr")
## write.csv(tt2[[1]], "Multinomial_Summary-Neg_ref_2yr.csv", row.names=FALSE)
## write.csv(tt2[[2]], "Multinomial_MC-Neg_ref_2yr.csv", row.names=FALSE)




# # This is in March after changing the strategy of not including the some predictors int he model.
# tt2 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv", 0.95, "FDR")
# write.csv(tt1[[1]], "Multinomial_Summary-NS_ref.csv", row.names=FALSE)
# write.csv(tt1[[2]], "Multinomial_MC-NS_ref.csv", row.names=FALSE)
 
# tt2 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv", 0.90, "FDR")
# tt3 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv", 0.85, "FDR")
# multinom_FDR = rbind(tt1, tt2, tt3)
# write.csv(multinom_FDR, "./species_data/ClimateNAData/Multinomial_summ_FDR_20180322.csv", row.names = FALSE)


#ss1 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv", 0.95, "non-FDR")
#ss2 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv", 0.90, "non-FDR")
#ss3 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data_SPSR_2.csv", 0.85, "non-FDR")
#multinom_NFDR = rbind(ss1, ss2, ss3)
#write.csv(multinom_NFDR, "./species_data/ClimateNAData/Multinomial_summ_NFDR_20180322.csv", row.names = FALSE)


### This is before March 2018

# tt1 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", 0.95, "FDR")
# tt2 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", 0.90, "FDR")
# tt3 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", 0.85, "FDR")


# ss1 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", 0.95, "non-FDR")

# ss2 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", 0.90, "non-FDR")

# ss3 = MultiMod("/srv/gspeedq32/body_size/species_data/ClimateNAData/meta_model_data.csv", 0.85, "non-FDR")




## Here pass the model and get the coefficient and 95% and 95% confidence interval 

MultiLogistic1 <- function(model, confintval, reference, VarName)
{
  ## Pass the model and generate the summary table here
  sum1 = summary(model)
  print(sum1)
  modDev = summary(model)$deviance
  modAIC = summary(model)$AIC
  
  ## Calculate the P value 
  model_z = (summary(model)$coefficients) / (summary(model)$standard.errors)
  
  model_p = (1 - pnorm(abs(model_z), 0, 1)) * 2
  
  ## Get the coefficients
  model_coef = exp(coef(model))
  ## Get the confidence intervals
  model_conf = exp(confint(model, level = confintval))
  
  ## confidence interval is displayed in a different format, so we need to transpose them
  dimnames( t(exp(confint(model, level = 0.95))[,,1]))
  ## initialise output matrix here
  TotRows = nrow(model_p)
  OpMat = matrix(0, nrow = TotRows, ncol = 23)
  OpMat[,1] = VarName
  OpMat[,2] = reference
  for (i in 1:TotRows)
  {
    Pval1 = model_p[i,]
    coef1 = model_coef[i,]
    ## confidence interval is a 3 dimentional array, 
    ## LowerConf1 = t(model_conf[,,1])[i,]
    ## UpperConf1 = t(model_conf[,,2])[i,]
    
    LowerConf1 = t(model_conf[,,i])[1,]
    UpperConf1 = t(model_conf[,,i])[2,]
    
    CompareTo = dimnames(model_coef)[[1]][i]
    DegreesofFreedom = nrow(model$fitted.values)
    OpMat[i,3:23] = c(CompareTo, coef1, LowerConf1, UpperConf1, Pval1, confintval, modDev, modAIC, DegreesofFreedom)
    
  }
  return(OpMat)
}

# op1 = MultiLogistic1(model1, 0.95, "NS", "bio1")

Multinom_plot <- function(prednames, model1, OpFileName, raw_dt)
{
  png(OpFileName, width = 10, height = 3, units = "in", res=300)
  #pdf(OpFileName, width = 10, height = 3)
  par(mfrow=c(1,3))
  #prednames = c('Mean_used', 'bio1_mean', 'bio1_breadth')
  for (count in 1:3) {
    # Rotate the variable names so we analyze each.
    if (count > 1) {
      prednames = c(prednames[2], prednames[3], prednames[1])
    }
    
    # Generate new data, holding two variables constant at their
    # means and covering the full range of the other.
    preddata = data.frame(
      col1=seq(from=min(raw_dt[[prednames[1]]]), to=max(raw_dt[[prednames[1]]]), length.out=1000),
      col2=rep(mean(raw_dt[[prednames[2]]]), 1000),
      col3=rep(mean(raw_dt[[prednames[3]]]), 1000)
    )
    colnames(preddata) = prednames
    head(preddata)
    
    # Calculate the multinomial parameter estimates for the new data.
    probs = predict(model1, newdata=preddata, 'probs')
    plot(preddata[[prednames[1]]], probs[,'NS'], ylim=c(0,1), type='l', xlab=prednames[1], ylab='Probability')
    lines(preddata[[prednames[1]]], probs[,'Neg'], col='red')
    lines(preddata[[prednames[1]]], probs[,'Pos'], col='blue')
  }
  
  dev.off()
  
}

