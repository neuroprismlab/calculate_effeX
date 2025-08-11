# Set parameters for calculate_effeX scripts

data_dir = 'tests/test_data'
intermediate_dir = 'tests/test_data/intermediates'
output_dir = 'tests/test_data/output'

final_output_file = 'braineffex_data'
final_output_path = file.path(output_dir, paste0(final_output_file, '_', Sys.Date(), '.RData'))

num_sdx_r2d = 2
alpha = 0.05

source('checker.R')
source('calc_d.R')
source('check_orientation.R')
source('clean_data.R')
source('helpers.R')

