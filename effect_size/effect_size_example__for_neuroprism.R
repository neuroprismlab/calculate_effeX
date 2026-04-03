# special script to run spatial extent results - 032526

library(devtools)
# install.packages("oro.nifti")
library(oro.nifti)
# install_github("neuroprismlab/BrainEffeX_utils")
library(BrainEffeX.utils)

# set params
percent <- '05' # USER-DEFINED

# set paths
data_dir <- paste0("/projects/neuroprism/effect_size/data/group_level/spatial_extent/",percent,"_percent/")
scripts_dir <- "/projects/neuroprism/effect_size/scripts/calculate_effeX/effect_size/"
intermediate_dir <- paste0("/projects/neuroprism/effect_size/data/combined_gl/intermediates/spatial_extent/",percent,"_percent/")
output_dir <- paste0("/projects/neuroprism/effect_size/data/combined_gl/output/spatial_extent/",percent,"_percent/")

# sub-paths
master_script <- paste0(scripts_dir,"master.R")
helper_dir <- paste0(scripts_dir,"scripts/")
template_filename <- paste0(scripts_dir,"data/template_nifti.nii.gz")

# make directories
dir.create(intermediate_dir)
dir.create(output_dir)

# Processing
source(master_script)
master(data_dir = data_dir,
       script_dir = helper_dir,
       intermediate_dir = intermediate_dir,
       output_dir = output_dir,
       template_filename = template_filename,
       final_output_file = paste0('braineffex_data_',percent,'_percent'))

# Visualization
# final_output_file <- file.path(output_dir, paste0('braineffex_data_',percent,'_percent', '_', Sys.Date(), '.RData'))
# load(final_output_file)
