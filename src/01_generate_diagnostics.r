# This file lists the species that we want to analyze. Contains 197 species.
sp_table_path = 'data_raw/species_table.csv'

# The directory with raw individual species data files.
data_dir = 'data_raw/bodymass_data/'

# Root directory for the modeling and diagnostics output.
output_dir = 'data_output/'

# Define the models that we want to fit.
model_specs = list(
  'm8' = massing ~ bio1 + bio12 + bio4 + bio15 + seasonality,
  'm10' = massing ~ bio1 + bio12 + bio15 + seasonality
)

# Define the models for which to do residual analyses.
res_models = list("m8", "m10")

# Make sure the output directories are ready.
clean_data_dir = paste0(output_dir, 'cleaned_bodymass_data/')
model_dir = paste0(output_dir, 'models/')
if (!dir.exists(output_dir)) { dir.create(output_dir) }
if (!dir.exists(clean_data_dir)) { dir.create(clean_data_dir) }
if (!dir.exists(model_dir)) { dir.create(model_dir) }

sp_table = read.csv(sp_table_path, header = TRUE, stringsAsFactors = FALSE)
if(!length(list.files(clean_data_dir))) 
  generateCleanedCSVFiles(data_dir, sp_table, clean_data_dir)

# results = Model_All(DataDir = clean_data_dir, sp_table, model_dir, model_specs, res_models)
