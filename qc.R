# Run QC on effect size FC matrices and output the plots to QC folder

fc_qc <- function(cleaned_data) {
    for (s in 1:dim(cleaned_data$study)[1]) {
        name <- cleaned_data$study$name[s]
        map_idx <- which(toupper(names(cleaned_data$effect_map)) == name)
        if ((cleaned_data$study$map_type[s] == "FC") & (cleaned_data$study$ref[s] == "Shen_268")) {
            plot_full_mat(cleaned_data$effect_map[[map_idx]]$orig_stat, "/work/neuroprism/effect_size/data/helper_data/map268_subnetwork.csv", export = TRUE, export_path = paste0("/work/neuroprism/effect_size/data/QC/", cleaned_data$study$name[s], "_orig_stat_mat.png"), show_plot = FALSE)
        }
        else if ((cleaned_data$study$map_type[s] == "FC") & (cleaned_data$study$ref[s] == "UKB_55")) {
            plot_full_mat(cleaned_data$effect_map[[map_idx]]$orig_stat, mapping_path = NA, export = TRUE, export_path = paste0("/work/neuroprism/effect_size/data/QC/", cleaned_data$study$name[s], "_orig_stat_mat.png"), show_plot = FALSE)
        }
    }
}


### QC function that creates HTML report (TESTING)
install.packages("rmarkdown")
install.packages("knitr")
install.packages("htmltools")
install.packages("base64enc")

library(rmarkdown)
library(knitr)
library(htmltools)
library(base64enc)

# Function to create and save HTML report
fc_qc <- function(cleaned_data) {
    # HTML content list
    html_content <- list()
    
    for (s in 1:nrow(cleaned_data$study)) {
        name <- cleaned_data$study$name[s]
        map_idx <- which(toupper(names(cleaned_data$effect_map)) == toupper(name))
        
        if ((cleaned_data$study$map_type[s] == "FC") & (cleaned_data$study$ref[s] %in% c("Shen_268", "UKB_55"))) {
            png_file <- tempfile(fileext = ".png")
            if (cleaned_data$study$ref[s] == "Shen_268") {
                plot_full_mat(cleaned_data$effect_map[[map_idx]]$orig_stat, 
                              "/work/neuroprism/effect_size/data/helper_data/map268_subnetwork.csv", 
                              export = TRUE, export_path = png_file, show_plot = FALSE)
            } else if (cleaned_data$study$ref[s] == "UKB_55") {
                plot_full_mat(cleaned_data$effect_map[[map_idx]]$orig_stat, 
                              mapping_path = NA, export = TRUE, export_path = png_file, show_plot = FALSE)
            }
            
            # Encode the image to base64
            encoded_image <- base64enc::dataURI(file = png_file, mime = "image/png")
            
            # Add the plot and its title to the HTML content list
            html_content <- c(html_content, list(h2(name), img(src = encoded_image, style = "width:25%;")))
        }
    }
    
    # Save the HTML content to a file
    html_file <- "/work/neuroprism/effect_size/data/QC/report.html"
    save_html(html_content, file = html_file)
    
    return(html_file)
}


