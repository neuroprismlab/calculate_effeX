# Sankey plot for EffeX studies

library(dplyr)
library(ggforce)
library(reshape2)

# data variable from combined_data should be called dataset
dataset <- data

# study dataframe renamed to data
data <- study

data$dataset <- toupper(data$dataset)

data$map_type <- toupper(data$map_type)

# add n manually for datasets
n_subs <- c("ABCD" = ">1500", "HBN" = ">400", "UKB" = ">35000", "HCP" = ">1000", "PNC" = ">400", "SLIM" = ">400")

# remove HCP_EP because small sample size
data <- data %>%
  filter(dataset != "HCP_EP")

data <- data %>%
  mutate(n_subs = n_subs[dataset])

# change hcp exp to hcp
# data <- data %>%
#   mutate(dataset = ifelse(dataset == "hcp_ep", "hcp", dataset))

data <- data %>%
  select(dataset, n_subs, map_type, orig_stat_type, category)

data <- data %>%
  group_by(dataset, n_subs, map_type, orig_stat_type, category) %>%
  summarise(n = n())

data <- reshape2::melt(data)
data <- gather_set_data(data, 1:5)


data <- data %>%
  mutate(dataset = reorder(dataset, -value))


specific_levels <- c("ABCD", "HCP", "PNC", "HBN", "SLIM", "UKB", ">1500", ">1000", ">400", ">35000","ACT","FC","t","t2","r")
# combine w remaining variables
combined_levels <- c(specific_levels, setdiff(levels(factor(data$y)), specific_levels))
data$y <- factor(data$y, levels = combined_levels)

# set colors
num_categories <- length(unique(data$y))
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_categories)








# Define the nodes that should be grey
grey_nodes <- c("FC", "ACT", "r", "t", "t2", "psychiatric", "cognitive", "biometric", "demographic", "socioeconomic")

# Create the plot
ggplot(data, aes(x, id = id, split = y, value = value)) +
  # The lines (parallel sets) will still be colored by dataset
  geom_parallel_sets(aes(fill = dataset), alpha = 0.3, axis.width = 0.1) +
  
  # The axes (nodes) will conditionally fill based on whether the node should be grey or not
  geom_parallel_sets_axes(aes(fill = ifelse(y %in% grey_nodes, "grey", y)), alpha = 0.7, axis.width = 0.1) +
  
  # Add the labels
  geom_parallel_sets_labels(colour = 'black', cex = 5, fontface = 'bold', angle = '0', hjust = 0, nudge_x = 0.08) +
  
  # Minimal theme and margins
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none", plot.margin = margin(0, 90, 0, -20)) +
  
  # Keep the lines (parallel sets) colored as per dataset
  scale_fill_manual(values = all_colors) +
  
  # Disable clipping so labels can extend beyond plot edges
  coord_cartesian(clip = "off")


# 
# # add corresponding n by dataset
# data2 <- data %>%
#   mutate(n_subs = n_subs[dataset])
# 
# # filter to only keep relevant columns
# data3 <- data2 %>%
#   select(dataset, n_subs, map_type, orig_stat_type)
#          
# data4 <- data3 %>%
#   group_by(dataset, map_type, orig_stat_type, n_subs) %>%
#   summarise(n = n())

# data5 <- reshape2::melt(data4)
# data6 <- gather_set_data(data5, 1:4)

ggplot(data6, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = dataset), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(aes(fill=y), alpha = 0.7, axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'black', cex = 5, fontface = 'bold', angle = '0', hjust = 0, nudge_x = 0.08) +
  theme_minimal() +
  theme(panel.grid = element_blank()) 





#########
data <- study
n_subs <- c("ABCD" = ">1500", "HBN" = ">400", "UKB" = ">35000", "HCP" = ">1000", "PNC" = ">400", "SLIM" = ">400")
data$dataset <- toupper(data$dataset)
data <- data %>%
  filter(dataset != "HCP_EP")
data <- data %>%
  mutate(n_subs = n_subs[dataset])

study_count <- data %>%
  count(dataset) %>%
  arrange(desc(n))

data$dataset <- factor(data$dataset, levels = study_count$dataset)
data$n_subs <- factor(data$n_subs, levels = unique(data$n_subs[order(match(data$dataset, levels(data$dataset)))]))

data <- data %>%
  select(dataset, n_subs, map_type, orig_stat_type, category) %>%
  group_by(dataset, n_subs, map_type, orig_stat_type, category) %>%
  summarise(n = n(), .groups = 'drop')

data <- reshape2::melt(data)

data$map_type <- as.factor(data$map_type)
data$category <- as.factor(data$category)


# Use gather_set_data from ggforce to prepare data for parallel sets plot
data <- gather_set_data(data, 1:5)

custom_colors <- c("ABCD" = "#66c2a5", "HBN" = "#fc8d62", "HCP" = "#8da0cb", "PNC" = "#e78ac3", "SLIM" = "#a6d854", "UKB" = "#ffd92f")
remaining_categories <- setdiff(unique(data$y), names(custom_colors))
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(remaining_categories))
remaining_colors <- setNames(color_palette, remaining_categories)
all_colors <- c(custom_colors, remaining_colors)

# Make factor and reorder in order to control plot order
# specific_levels <- c("HBN", "PNC", "SLIM", "HCP", "ABCD", "UKB", ">400", ">1000", ">1500", ">35000", "ACT", "FC", "t", "t2", "r")
# combined_levels <- c(levels(data$dataset), ">400", ">1000", ">1500", ">35000", "ACT", "FC", "t", "t2", "r")
# data$y <- factor(data$y, levels = combined_levels)

ggplot(data, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = dataset), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(aes(fill = y), alpha = 0.7, axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'black', cex = 5, fontface = 'bold', angle = '0', hjust = 0, nudge_x = 0.08) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none", plot.margin = margin(0, 90, 0, -20)) +
  scale_fill_manual(values = all_colors) +
  coord_cartesian(clip = "off")



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
  scale_fill_manual(values = color_palette) +
  coord_cartesian(clip = "off")
