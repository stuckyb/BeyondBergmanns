# This file lists the species that we want to analyze. Contains 197 species.
sp_table_path = here('data_raw/species_table.csv')

# The directory with raw individual species data files.
data_dir = here('data_raw/bodymass_data/')

# Root directory for the modeling and diagnostics output.
output_dir = here('data_output/')


# Define model specification with sex as covariate for 2 year 
model_specs = list(
  'm8' = massing ~ bio1_2yr + bio12_2yr + bio4_2yr + bio15_2yr + seasonality + sex_cor,
  'm10' = massing ~ bio1_2yr + bio12_2yr + bio15_2yr + seasonality + sex_cor
)


# Define the models for which to do residual analyses.
res_models = list("m8", "m10")

# Make sure the output directories are ready.
# clean_data_dir = here('data_output/cleaned_bodymass_data/')
# This directory where sex and lifestage columns are added in the species list. 
clean_data_dir = here('data_output/cleaned_bodymass_data/')

## For Linux
# model_dir = here('data_output/models/')

## For windows
model_dir = paste(here('data_output/models/'), "/", sep = "")


if (!dir.exists(output_dir)) { dir.create(output_dir) }
if (!dir.exists(clean_data_dir)) { dir.create(clean_data_dir) }
if (!dir.exists(model_dir)) { dir.create(model_dir) }

sp_table = read.csv(sp_table_path, header = TRUE, stringsAsFactors = FALSE)
if(!length(list.files(clean_data_dir))) 
  generateCleanedCSVFiles(data_dir, sp_table, clean_data_dir)

## Run the models with bioclimatic predictors, season of collection and sex as predictor
## results = Model_All(DataDir = clean_data_dir, sp_table, model_dir, model_specs, res_models, "all")

