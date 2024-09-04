# Set parameters for calculate_effeX scripts

data_dir = '/work/neuroprism/effect_size/data/group_level'
script_dir = '/home/h.shearer/hallee/calculate_effeX/combine_gl'
intermediate_dir = '/work/neuroprism/effect_size/data/combined_gl/intermediates'
output_dir = '/work/neuroprism/effect_size/data/combined_gl/output'

final_output_file = 'combined_data'
final_output_path = file.path(output_dir, paste0(final_output_file, '_', Sys.Date(), '.RData'))

num_sdx_r2d = 2
alpha = 0.05