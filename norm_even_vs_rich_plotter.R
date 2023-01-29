###################################################
############ Austin Marshall and ChatGPT ##########
###################################################

# Import necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)

# read in tsv count file
counts_data <- read.delim("emu-combined-abundance-species-counts.tsv", sep = "\t", header = TRUE, row.names = 1) %>%
  mutate_all(~replace(., is.na(.), 0))

# calculate species richness
species_richness <- apply(counts_data, 2, function(x) sum(x > 0))

# calculate normalized mean evenness
normalized_mean_evenness <- apply(counts_data, 2, function(x) diversity(x, "shannon")/log(sum(x > 0)))

# create data frame for plotting
plot_data <- data.frame(sample = colnames(counts_data),
                        species_richness = species_richness,
                        normalized_mean_evenness = normalized_mean_evenness)

# create plot
ggplot(plot_data, aes(x = species_richness, y = normalized_mean_evenness, color = sample)) +
  geom_point() +
  xlab("Species Richness") +
  ylab("Normalized Mean Evenness") +
  ggtitle("16S Aerosol Sequencing with Varying Input Conc. and Cycling")
