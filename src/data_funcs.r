
source('src/collin_diags.r')


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
        CurRow = c(SpName)

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

