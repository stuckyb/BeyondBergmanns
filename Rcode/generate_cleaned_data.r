
# This file lists the species that we want to analyze. Contains 197 species
species_path = "data_raw/species_table.csv"

# The directory with individual species data files.
DataDir = "data_raw/bodymass_data/"



#------
# Generates cleaned body mass data files.  Filters out records with spurious
# temperature and/or body size data.
#
# DataDir (str): The directory with individual species data files.
# RangeFileName (str): File that lists the species that we want to analyze.
# output_dir: Where to save the cleaned CSV files.
# FromRecord (integer): Start of the row range to analyze.
# ToRecord (integer): End of the row range to analyze (inclusive).
#------
generateCleanedCSVFiles = function(DataDir, species_path, output_dir, FromRecord, ToRecord)
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


n1 = generateCleanedCSVFiles(DataDir, species_path, 'data_output/cleaned_bodymass_data/', 1, 198)

