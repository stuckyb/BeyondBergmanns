library(tidyverse)
#install.packages("compositions")
library(compositions)
library(vegan)


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