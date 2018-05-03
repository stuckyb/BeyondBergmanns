#setwd("C:/Narayani/Projects/BodySize_BergmansRule/test_data/")


# Load the conditioning/collinearity diagnostics functions.  (You might need to
# change the path, depending on where things are located on your system.)
source('../conditioning_diagnostics/conditioning_diagnostics.r')


## BuildModels for each species
## take the required data from the complete data for building models
## Required columns are ("month", "decimallongitude", "decimallatitude", "lifestage", "year", "scientificname", "lengthunitsinferred", "massing", "lengthinmm", "sex", "bio1", "bio12")

#setwd("/srv/gspeedq32/narayani/body_size/test_data/")


#SpListFile = "SpList.csv"
#RangeFileName = "SpList_massRange.csv"

# This file lists the species that we want to analyze. Contains 197 species
RangeFileName = "./species_data/SpList_sorted185_n1.csv"

## This file contains only 187 species.
#RangeFileName = "./species_data/SpList_sorted185_n120170706_c.csv" 

# The directory with individual species data files.
DataDir = "./species_data/ClimateNAData/output/"
# Directory for the modeling results.
#resDir = "./species_data/ClimateNAData/results/"
#resDir = "./species_data/ClimateNAData/results_new_0.15_2.25/"

## This folder contains the results with standardised predictors and body mass as response with actual values
#resDir = "./species_data/ClimateNAData/results_new_0.25_2_trim_mean/"

## This forlder contains the results with centered predictors and standardised body mass as response
#resDir = "./species_data/ClimateNAData/results_new_0.25_2_trim_mean_center_pred_std_response/"

## This folder contains the results with no transformation of predictors and standardised body mass as response
## resDir = "./species_data/ClimateNAData/results_new_0.25_2_trim_mean_act_pred_std_response/"

## This folder contains the results with for confirming the initial results, trimmed mean 0.25_2, raw predictors and standardised response, Date 20/2/2018
resDir = "./species_data/ClimateNAData/results_valid_m8_m10_models_SPSR/"


# # Define the models that we want to fit.
model_specs = list(
  #'m1' = massing ~ bio1 + bio12,
  #'m2' = massing ~ bio1 + bio12 + bio1:bio12,
  #'m3' = massing ~ bio1 + bio12 + seasonality,
  #'m4' = massing ~ bio1 + bio12 + seasonality + bio1:bio12,
  #'m5' = massing ~ bio1 + bio12 + bio4 + bio15 + bio1:bio12,
  # 'm6' = massing ~ bio1 + bio12 + bio4 + bio15 + seasonality + bio1:bio12,
  #'m7' = massing ~ bio1 + bio12 + bio15,
  'm8' = massing ~ bio1 + bio12 + bio4 + bio15 + seasonality,
  #'m9' = massing ~ bio1 + bio12 + bio4 + bio15 + seasonality + bio1:bio12 + bio1:bio4 + bio12:bio15 + bio4:bio15,
  'm10' = massing ~ bio1 + bio12 + bio15 + seasonality
  #'m11' = massing ~ bio1 + bio12 + bio4 + seasonality,
  #'m12' = massing ~ bio1 + bio12 + bio4 + bio15
)



# Define the models that we want to fit. Want to fit only m8, decided that m8 is the model which fits the data well. 
# model_specs = list(
#  'm8' = massing ~ bio1 + bio12 + bio4 + bio15 + seasonality
#  )


# Define the models for which to do residual analyses.
res_models = list("m8", "m10")


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
# Standardizes one or more columns in a dataframe.
# predictors are centered and not scaled. 
#----

standardize_cols = function(df, colnames) {
  for (colname in colnames) {
    df[[colname]] = as.vector(scale(df[[colname]], scale=TRUE, center=TRUE))
  }
  
  return(df)
}
# Simple test of standardize_cols().
(testdata = data.frame(
  a=rnorm(10, mean=10, sd=5), b=rnorm(10, mean=50, sd=10), c=rnorm(10, mean=100, sd=20)
))
colMeans(testdata)
sd(testdata$a)
sd(testdata$b)
sd(testdata$c)
testdata = standardize_cols(testdata, c('a', 'c'))
colMeans(testdata)
sd(testdata$a)
sd(testdata$b)
sd(testdata$c)



PlotResiduals <- function(CurModel, SpName, mod_name, resDir)
{
#	jpegfileName = paste("./species_data/ClimateNAData/results/", SpName, "_", mod_name, "_res.jpg", sep = "")
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
#	resFileName = paste("./species_data/ClimateNAData/results/", SpName, "_", mod_name, "_res.csv", sep = "")
	resFileName = paste(resDir, SpName, "_", mod_name, "_res.csv", sep = "")	
	print(resFileName)
	names(mat1)[1:2] = c("residuals", "fitted.values")
	write.csv(mat1, resFileName, row.names=FALSE)
	
}



#------
# Runs model fits for a range of records from the master species list file.
#
# DataDir (str): The directory with individual species data files.
# RangeFileName (str): File that lists the species that we want to analyze.
# resDir (str): Directory for the modeling results.
# FromRecord (integer): Start of the row range to analyze.
# ToRecord (integer): End of the row range to analyze (inclusive).
#------
Model_All <- function(DataDir, RangeFileName, resDir, FromRecord, ToRecord)
{
#	t1 = read.csv(SpListFile, header=TRUE)

	RangeTbl = read.csv(RangeFileName, header=TRUE, stringsAsFactors = FALSE)
	#t1 = RangeTbl[,1]
	TotSps = dim(RangeTbl)[1]
	num_rescols = 1 + length(model_specs) * 7
	ResMat = matrix(0, nrow=0, ncol = num_rescols)
	MeanTbl = matrix(0,nrow=TotSps, ncol = 9)
	CollinTbl = matrix("", nrow=0, ncol=13)
	preds_list = c(
	    '(Intercept)', 'bio1', 'bio12', 'bio4', 'bio15', 'seasonalitySpring',
		'seasonalitySummer', 'seasonalityWinter'
	)
	#for (i in 1:TotSps)
	for (i in FromRecord:ToRecord)
	{
		SpName = trimws(as.character(RangeTbl[i,1]), which= "both")
		InpFileName = paste(DataDir, "Env_", SpName, ".csv", sep = "")
		
		print(paste(SpName, InpFileName, sep = ":"))
		
		Tbl = read.csv(InpFileName, header=TRUE, stringsAsFactors=FALSE, fileEncoding="UTF-8" )
		print(dim(Tbl))

		# Keep only the columns with data that we (might) want to use
		# in the models:
		# lengthunitsinferred (bool): ??
		# massing: mass in grams
		# bio1: Annual mean temp.
		# bio12: Annual total precip.
		# bio4: Annual temperature seasonality.
		# bio15: Annual precipitation seasonality.
		# seasonality: Season (categorical).
		#Tbl2 = Tbl[, c("month", "decimallongitude", "decimallatitude", "lifestage", "year", "scientificname", "lengthunitsinferred", "massing", "lengthinmm", "sex", "bio1", "bio12", "bio4", "bio15", "seasonality")]
		Tbl2 = Tbl
		
		
		# These define the mass range for records that we want to include in the model.
		MinMass = as.numeric(RangeTbl[i, 5])
		MaxMass = as.numeric(RangeTbl[i, 6])
		
		## This is Literature mean. 
		MidPoint = MinMass + (MaxMass - MinMass) / 2
		
		## Remove data rows with temperature is too low. Specially to remove data with -9999 values in bio1
		NoDataIndices = which(Tbl2[['bio1']] <= -50)
		print(paste("Total -9999 rows ", length(NoDataIndices)))
		if (length(NoDataIndices) > 0 )
		{
			Tbl2 = Tbl2[-NoDataIndices,]
		}
		 print(dim(Tbl2))
		 
		 
		## Find the truncated mean 10% of the data from upper and lower range removed. 
		## Get the lower and upper bounds from the data, so that SD can be calculated. 
		Tbl2 = Tbl2[order(Tbl2$massing), ]
		trim_pct = 0.1
		TrimmedMean = mean(Tbl2$massing, na.rm = TRUE, trim = trim_pct)
		
		MeanTbl[i,1:7] = c(SpName, Tbl2[round(nrow(Tbl2)*trim_pct), c('massing')], Tbl2[round(nrow(Tbl2)*(1-trim_pct)), c('massing')], TrimmedMean, MinMass, MaxMass, MidPoint)
		
		## Changing the Midpoint to trimmed mean
		MidPoint = TrimmedMean		
		MinMass = MidPoint * 0.25
		MaxMass = MidPoint * 2.0
		
		ValidIndices = which(Tbl2[['massing']] >= MinMass & Tbl2[['massing']] <= MaxMass)
		Tbl3 = Tbl2[ValidIndices,]
		Tbl3 = Tbl3[Tbl3$seasonality != '',]

		# Initialize the output results row.
		CurRow = c(SpName)

		MeanTbl[i,8:9] = c(min(Tbl3$massing, na.rm = TRUE), max(Tbl3$massing, na.rm = TRUE))
		
		# Do each model fit and get the fit statistics. Do not standardize_cols for APSR, center columns for CPSR. 
		Tbl3_stdzd = standardize_cols(Tbl3, c('bio1', 'bio12', 'bio4', 'bio15'))
		Tbl3_stdzd$massing = scale(Tbl3_stdzd$massing)

		for (mod_name in names(model_specs)) {
		  model_spec = model_specs[[mod_name]]
		  
		  collin_res = collin_diags(model_spec, Tbl3)
		  #print(str(collin_res))
		  c_num = collin_res$c_num
		  #c_num = -9999
		  #print(paste('CONDNUM:', c_num))
		  
		  ncols = ncol(collin_res$VDPs)
		  newrow = rep(0, length(preds_list))
		  row_exist = which(preds_list %in% row.names(collin_res$VDPs))
		  newrow[row_exist] = collin_res$VDPs[,ncols]
		  CollinTbl = rbind(
			CollinTbl,
			c(SpName, mod_name, as.character(collin_res$c_indices[1:3]), as.character(newrow))
		  )

		  if (c_num > 30) {
		    fname = paste('collin_results/', SpName, '-', mod_name, '.txt', sep='')
		    write_collin_report(fname, collin_res)
		  }

		  # Standardardize the continuous variables before model fitting.

		  mod = lm(model_spec, data=Tbl3_stdzd)

		  ## This is to plot the residuals and save the residuals and fitted values in a text file. 
		  if (length(which(mod_name %in% res_models)) > 0 )
		  {			
			PlotResiduals(mod, SpName, mod_name, resDir)
		  }

		  mod_summary = GetModelSummary (mod, resDir, SpName, mod_name)
		  
		  CurRow = c(CurRow, mod_summary, c_num, AIC(mod))
		  print(CurRow)

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
	## Centered Predictors Standardized Response CPSR
	#write.csv(res_df, paste(resDir, '../model_summaries-climateNA-0.25_2_trim_mean_CPSR.csv', sep=''), row.names=FALSE)
	#write.csv(MeanTbl, paste(resDir, '../MeanTable_CPSR.csv', sep=''), row.names=FALSE) 
	## Actual Predictors Standardized Response APSR
	write.csv(res_df, paste(resDir, '/model_summaries-climateNA-0.25_2_trim_mean_m8_m10_SPSR.csv', sep=''), row.names=FALSE)
	write.csv(MeanTbl, paste(resDir, '/MeanTable_0.25_2_trim_mean_m8_m10_SPSR.csv', sep=''), row.names=FALSE) 
	CollinTbl = as.data.frame(CollinTbl)
	names(CollinTbl) = c('SpeciesName', 'mod', 'ci_1', 'ci_2', 'ci_3', preds_list)
	write.csv(CollinTbl, paste('collin_results/', '0_2models_condition_m8_m10_summaries.csv', sep=''), row.names=FALSE) 
	
	
	#write.csv(res_df, paste(resDir, '../model_summaries-climateNA_tt.csv', sep=''), row.names=FALSE)
	return(ResMat)
}

#n1 = Model_All()

GetModelSummary <- function(model, resDir, SpName, modelName)
{
	## Get the coefficients and store it in matrix
	OpMat1 = summary(model)$coefficients
	
	## Get the string of significant variables in variable SigVar
	#j = which(OpMat1[-1,4] < 0.05)
	#SigVarName = paste(names(j), sep= "", collapse = " ; ")
	## To add (+) or (-) sign in the significant variable. 
	pos = which(OpMat1[-1,4] < 0.05 & OpMat1[-1,1] > 0)
	neg = which(OpMat1[-1,4] < 0.05 & OpMat1[-1,1] < 0)
	SigVarNamePos = ''
	SigVarNameNeg = ''
	SigVarName = ''
	if (length(pos) > 0)
	{
		SigVarNamePos = paste(names(pos), rawToChar(as.raw(40)), rawToChar(as.raw(43)), rawToChar(as.raw(41)), collapse = "; ", sep = "")
		SigVarName = SigVarNamePos
		#print("In positive")		
	} 
	
	#print(SigVarNamePos)
	if (length(neg) > 0)
	{
		SigVarNameNeg = paste(names(neg), rawToChar(as.raw(40)), rawToChar(as.raw(45)), rawToChar(as.raw(41)), collapse = "; ", sep = "")
		SigVarName = SigVarNameNeg
		#print("In negative")
	} 
	#print(SigVarNameNeg)
	
	if (nchar(SigVarNamePos) & nchar(SigVarNameNeg))
	{
		#print("Here")
		SigVarName = paste(SigVarNamePos, SigVarNameNeg, sep = "; ")
	}
	
	
	#print(SigVarName)	
	#print(length(SigVarName))
	
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
	
	#print("Before writing")
	## save the models in results directory
	OpFileName = paste(resDir, SpName, "_", modelName, ".csv", sep = "")
	write.csv(OpMat1, OpFileName, row.names=FALSE)	
	
	## Now get the last four rows from each model and transpose it to columns and attach it to output matrix
	TotRows= dim(OpMat1)[1]
	nt1 = OpMat1[(TotRows - 3):TotRows, 2]	
	#print(nt1)
	nt1 = c(nt1, SigVarName)
	#print(nt1)
	#return(OpMat1)
	return(nt1)

}

# test1 = GetModelSummary(m2, resDir, SpName, modelName)

n1 = Model_All(DataDir, RangeFileName, resDir, 1, 198)
#n1 = Model_All(DataDir, RangeFileName, resDir, 1, 5)
#n1 = Model_All(DataDir, RangeFileName, resDir, 105, 111)
#n1 = Model_All(DataDir, RangeFileName, resDir, 175, 175)


