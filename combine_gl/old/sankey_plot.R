# Sankey plot for EffeX studies

library(dplyr)
library(ggforce)
library(reshape2)
library(RColorBrewer)

load("/work/neuroprism/effect_size/effect_maps_clean.RData")

data <- study

## Post-hoc meta-data corrections

# add n manually for datasets
n_subs <- c("ABCD" = ">1500", "HBN" = ">400", "UKB" = ">35000", "HCP" = ">1000", "PNC" = ">400", "SLIM" = ">400")

# remove IMAGEN
data <- data %>%
  filter(dataset != "IMAGEN")

# rename codes
data$code[data$code=='CBCL' | data$code=='CBCL (follow-up)'] <-'psychiatric'
# data$code[data$code=='CBCL' | data$code=='CBCL (follow-up)'] <-'psychiatric (CBCL)'
data$code[data$code=='clinical'] <- 'psychiatric'
data$code[data$code=='cognitive (IQ)'] <- 'cognitive'
data$code[data$code=='demographic (age)'] <- 'demographic'
data$code[data$code=='demographic (gender)'] <- 'demographic'
data$orig_stat_type[data$orig_stat_type=='d'] <- 't'
# data$map_type[data$map_type=='ACT'] <- 'Activation'

# change hcp exp to hcp
data <- data %>%
  mutate(dataset = ifelse(dataset == "HCPXEP", "HCP", dataset))

# add corresponding n by dataset
data <- data %>%
  mutate(n_subs = n_subs[dataset])

# filter to only keep relevant columns
data <- data %>%
  select(dataset, n_subs, map_type, orig_stat_type, code)

data <- data %>%
  group_by(dataset, n_subs, map_type, orig_stat_type, code) %>%
  summarise(n = n())

# data <- data %>%
#   group_by(dataset, map_type, orig_stat_type, n_subs, code) %>%
#   summarise(n = n())

data <- reshape2::melt(data)
data <- gather_set_data(data, 1:5)

# custom_colors <- c("ABCD" = "#66c2a5", "HBN" = "#fc8d62", "HCP" = "#8da0cb", "PNC" = "#e78ac3", "SLIM" = "#a6d854", "UKB" = "#ffd92f")
# remaining_categories <- setdiff(unique(data$y), names(custom_colors))
# color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(remaining_categories))
# remaining_colors <- setNames(color_palette, remaining_categories)
# all_colors <- c(custom_colors, remaining_colors)

# make factor and reorder in order to control plot order
specific_levels <- c("HBN", "PNC", "SLIM", "HCP", "ABCD", "UKB", ">400", ">1000", ">1500", ">35000","ACT","FC","t","t2","r")
# combine w remaining variables
combined_levels <- c(specific_levels, setdiff(levels(factor(data$y)), specific_levels))
data$y <- factor(data$y, levels = combined_levels)

# set colors
num_categories <- length(unique(data$y))
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_categories)

ggplot(data, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = dataset), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(aes(fill=y), alpha = 0.7, axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'black', cex = 5, fontface = 'bold', angle = '0', hjust = 0, nudge_x = 0.08) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none", plot.margin = margin(0, 90, 0, -20)) +
  scale_fill_manual(values = all_colors) +
  coord_cartesian(clip = "off")
