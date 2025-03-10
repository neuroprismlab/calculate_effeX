library(ggforce)
# ?geom_parallel_sets

# from example: keeping link color constant https://stackoverflow.com/questions/57506648/is-there-an-option-to-keep-the-same-color-for-nodes-with-the-same-label-sankey
# data <- reshape2::melt(Titanic)
# data <- gather_set_data(data, 1:4)
# 
# ggplot(data, aes(x, id = id, split = y, value = value)) +
#   geom_parallel_sets(aes(fill = Sex), alpha = 0.3, axis.width = 0.1) +
#   geom_parallel_sets_axes(aes(fill=y), axis.width = 0.1) +
#   geom_parallel_sets_labels(colour = 'white')


# create simple dataset that replicates titanic dataset

library(dplyr)
data2 <- data.frame(
  dataset = c("ABCD", "ABCD", "ABCD", "hbn", "hbn", "1. hcp", "1. hcp", "pnc", "pnc", "pnc"),
  type = c("1. fc", "1. fc", "1. fc", "1. fc", "1. fc", "1. fc", "act", "1. fc",  "1. fc", "1. fc"),
  statistic = c("r", "r", "t", "r", "r", "r", "t2", "r",  "r", "r"),
  condition1 = c("rest", "rest", "rest", "rest", "rest", "rest", "rest", "rest", "rest", "rest"),
  condition2 = c("age", "bmi", "iq","nihtoolbox", "srs", "age","age", "age", "age", "age"),
  n_subs = c('>1500', '>1500', '>1500', '>400', '>400', '>1000', '>1000', '>400', '>400', '>400'),
  condition_category = c("1. demographic", "physical", "cognitive","cognitive", "psychiatric", "1. demographic","1. demographic", "1. demographic", "1. demographic", "1. demographic")
  )

# TODO: the numbers are a workaround to plit in the order I want, since geom_parallel_sets plots in alphabetical order
# TODO: figure out reordering for plot. The below releveling does not work...
# data2$dataset <- factor(data2$dataset, levels = c("hcp", "ABCD", "hbn", "pnc"))


# for data2, tabulate number of studies in each combination of categories then make it long - needed for sankey plotting
data2 <- data2 %>% 
  group_by(dataset, n_subs, type, statistic, condition_category) %>% 
  summarise(n = n())
data2 <- reshape2::melt(data2)
data2 <- gather_set_data(data2, 1:5)

ggplot(data2, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = dataset), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(aes(fill=y), axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'white')






# # OLD EXAMPLE:
#
# library(networkD3)
# 
# # Create a data frame of nodes
# nodes <- data.frame(
#   name = c("45 Studies", "Dataset 1", "Dataset 2", "Age Range 1", "Age Range 2")
# )
# 
# # Create a data frame of links
# links <- as.data.frame(matrix(c(
#   0, 1, 22, 1, # 22 studies from "45 Studies" to "Dataset 1", group 1
#   0, 2, 23, 2, # 23 studies from "45 Studies" to "Dataset 2", group 2
#   1, 3, 15, 1, # 15 studies from "Dataset 1" to "Age Range 1", group 1
#   1, 4, 7, 1,  # 7 studies from "Dataset 1" to "Age Range 2", group 1
#   2, 3, 10, 2, # 10 studies from "Dataset 2" to "Age Range 1", group 2
#   2, 4, 13, 2  # 13 studies from "Dataset 2" to "Age Range 2", group 2
# ), ncol = 4, byrow = TRUE))
# 
# # Name the columns in the links data frame
# names(links) <- c("source", "target", "value", "group")
# 
# # Create a color scale
# colourScale <- brewer.pal(2, "Set1")
# 
# # Create the Sankey plot
# sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
#               Target = "target", Value = "value", NodeID = "name",
#               units = "Studies", fontSize = 12, nodeWidth = 30)
# 
# sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
#               Target = "target", Value = "value", NodeID = "name",
#               units = "Studies", fontSize = 12, nodeWidth = 30,
#               LinkGroup = "group")
# 

