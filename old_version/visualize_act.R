# save mean act
library(oro.nifti)

act_dir<-"/Users/stephanienoble/Library/CloudStorage/GoogleDrive-stephanie.noble@yale.edu/My Drive/Lab/More/Tasks-Ongoing/K99/Effect_Size/data/datasets/misc/hcp maps/HCP/hcp task activ"

# list all files that end with ".nii.gz" in folder "/Users/stephanienoble/Library/CloudStorage/GoogleDrive-stephanie.noble@yale.edu/My Drive/Lab/More/Tasks-Ongoing/K99/Effect_Size/data/datasets/misc/hcp maps/HCP"

# act_studies <- names(d_master__reordered)[grep("act", names(d_master__reordered))]

act_files<-list.files(act_dir, pattern = ".nii.gz", full.names = TRUE)

mat_to_plot<-vector(mode = "list", length=length(act_files))
names(mat_to_plot)<-act_files
# read files into mat_to_plot in a loop
for (i in 1:length(act_files)) {
# for (this_file in list.files(act_dir, pattern = ".nii.gz", full.names = TRUE)) {
    mat_to_plot[[act_files[i]]] <- readNIfTI(act_files[i])
}

# average across list of matrices in mat_to_plot
mat_to_plot__mean <- Reduce("+", mat_to_plot)/length(mat_to_plot)

writeNIfTI(mat_to_plot__mean, paste0(act_dir, "/mean_act")) # automatically adds .nii.gz

# Visualize - using BISWeb

# image(mat_to_plot__mean, useRaster = TRUE)
# plot nifti image
ortho2(mat_to_plot__mean)