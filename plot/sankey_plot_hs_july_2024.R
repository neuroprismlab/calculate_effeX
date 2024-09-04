# Sankey plot for EffeX studies

library(dplyr)
library(ggforce)
library(reshape2)

data <- study

# add n manually for datasets
n_subs <- c("ABCD" = ">1500", "HBN" = ">400", "UKB" = ">35000", "HCP" = ">1000", "PNC" = ">400", "SLIM" = ">400")

# remove IMAGEN
data <- data %>%
  filter(dataset != "IMAGEN")

# change hcp exp to hcp
data <- data %>%
  mutate(dataset = ifelse(dataset == "HCPXEP", "HCP", dataset))


# add corresponding n by dataset
data2 <- data %>%
  mutate(n_subs = n_subs[dataset])

# filter to only keep relevant columns
data3 <- data2 %>%
  select(dataset, n_subs, map_type, orig_stat_type)
         
data4 <- data3 %>%
  group_by(dataset, map_type, orig_stat_type, n_subs) %>%
  summarise(n = n())

data5 <- reshape2::melt(data4)
data6 <- gather_set_data(data5, 1:4)

ggplot(data6, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = dataset), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(aes(fill=y), alpha = 0.7, axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'black', cex = 5, fontface = 'bold', angle = '0', hjust = 0, nudge_x = 0.08) +
  theme_minimal() +
  theme(panel.grid = element_blank()) 



