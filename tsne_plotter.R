library(vegan)
library(Rtsne)

a = read.delim(file = "emu-combined-abundance-species-counts.tsv", check.names = FALSE) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(species, `8_9_2005`, `10_11_2009`, `6_2_2010_A`, `6_2_2010_B`, `10_25_2014`, `5_11_2016`, `8_26_2016_A`,
         `8_26_2016_B`, `3_19_2017`, `7_3_2017`, `9_17_2017`, `1_13_2018`, `8_3_2018`, `1_14_2019_A`, `1_14_2019_B`,
         `6_3_2019`, `1_7_2020_A`, `1_7_2020_B`, `4_7_2020_A`, `4_7_2020_B`, `8_1_2020_A`, `8_1_2020_B`, `1_14_2021_A`,
         `1_14_2021_B`, `7_10_2021_A`, `7_10_2021_B`, `7_10_2021_C`, `7_10_2021_D`, `10_1_2021`, `1_24_2022_A`, `1_24_2022_B`, 
         `8_20_2022_A`, `8_20_2022_B`, `8_20_2022_C`, `1_11_2023_A`, `1_11_2023_B`, `1_11_2023_C`, `1_11_2023_D`, `1_11_2023_E`)

b <- a %>%
  filter(!grepl("\\[", species))

colnames(b)[1] = "rowname"

otu_data = b %>%
  column_to_rownames("rowname") %>%
  t()

# Convert the OTU data to numeric (in case there are non-numeric values)
otu_data <- as.matrix(otu_data)
otu_data <- apply(otu_data, 2, as.numeric)

# Handle missing values (if any)
# If you want to remove samples with missing values, use:
# otu_data <- na.omit(otu_data)
# If you want to impute missing values, use an appropriate method (e.g., mean imputation, k-nearest neighbors, etc.)

# Cumulative Sum Scaling (CSS) filtering using vegan package
filtered_data <- decostand(otu_data, method = "normalize")

# Remove duplicate rows (if any)
filtered_data <- filtered_data[!duplicated(filtered_data), ]

# t-SNE analysis with adjusted perplexity
# Replace the "your_perplexity_value" with an appropriate value (e.g., between 5 and nrow(filtered_data) - 1)
perplexity_value <- 10
set.seed(123)
tsne_result <- Rtsne(filtered_data, perplexity = perplexity_value)

# Create a t-SNE plot
# Replace "your_labels.csv" with the path to a file containing sample labels (one per line)
sample_labels <- read.csv("same2022_sample2plot.csv", header = TRUE)
sample_labels <- as.character(sample_labels$Plot_name)

# Plot t-SNE
plot(tsne_result$Y, col = "blue", pch = 19, xlab = "t-SNE 1", ylab = "t-SNE 2", cex = 2.25, main = "16S SLR Microbiomes")

# Add sample labels to the plot
text(tsne_result$Y, labels = sample_labels, col = "black", cex = 1.25, pos = 3, )

# Save the plot (optional)
# Replace "tsne_plot.png" with the desired output file name and format (e.g., "tsne_plot.pdf")
# savePlot("tsne_plot.png")

a = png("dust_16S_tsne_plot.png", width = 8, height = 6, units = "in", res = 600)
dev.off()

###############___NMDS___##################################################################

# Load required libraries
install.packages("vegan")  # If you haven't installed it before
library(vegan)

a = read.delim(file = "emu-combined-abundance-species-counts.tsv", check.names = FALSE) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(species, `8_9_2005`, `10_11_2009`, `6_2_2010_A`, `6_2_2010_B`, `10_25_2014`, `5_11_2016`, `8_26_2016_A`,
         `8_26_2016_B`, `3_19_2017`, `7_3_2017`, `9_17_2017`, `1_13_2018`, `8_3_2018`, `1_14_2019_A`, `1_14_2019_B`,
         `6_3_2019`, `1_7_2020_A`, `1_7_2020_B`, `4_7_2020_A`, `4_7_2020_B`, `8_1_2020_A`, `8_1_2020_B`, `1_14_2021_A`,
         `1_14_2021_B`, `7_10_2021_A`, `7_10_2021_B`, `7_10_2021_C`, `7_10_2021_D`, `10_1_2021`, `1_24_2022_A`, `1_24_2022_B`, 
         `8_20_2022_A`, `8_20_2022_B`, `8_20_2022_C`, `1_11_2023_A`, `1_11_2023_B`, `1_11_2023_C`, `1_11_2023_D`, `1_11_2023_E`)

b <- a %>%
  filter(!grepl("\\[", species))

colnames(b)[1] = "rowname"

otu_data = b %>%
  column_to_rownames("rowname") %>%
  t()

# Convert the OTU data to numeric (in case there are non-numeric values)
otu_data <- as.matrix(otu_data)
otu_data <- apply(otu_data, 2, as.numeric)

# Cumulative Sum Scaling (CSS) filtering using vegan package
filtered_data <- decostand(otu_data, method = "normalize")

# Remove duplicate rows (if any)
filtered_data <- filtered_data[!duplicated(filtered_data), ]

# Handle missing values (if any)
# If you want to remove samples with missing values, use:
# otu_data <- na.omit(otu_data)
# If you want to impute missing values, use an appropriate method (e.g., mean imputation, k-nearest neighbors, etc.)

# NMDS analysis
nmds_result <- metaMDS(filtered_data, distance = "bray")

# Replace "your_labels.csv" with the path to a file containing sample labels (one per line)
sample_labels <- read.csv("ISS_OTU_species.csv", header = TRUE)
sample_labels <- as.character(sample_labels$Plot_name)

# Create an NMDS plot
plot(nmds_result, type = "n")  # "n" means no points initially

# Add points to the plot
points(nmds_result, col = "blue", pch = 19)

# Add sample labels to the plot
text(nmds_result, labels = sample_labels, col = "red", cex = 0.7, pos = 3)

# Add a legend to the plot (if needed)
legend("topright", legend = unique(sample_labels), col = unique(sample_labels), pch = 19)

# Add a title to the plot
title("NMDS Plot of OTU Data")

# Save the plot at a higher resolution (e.g., 300 dpi) as a PNG image
# Replace "nmds_plot.png" with the desired output file name and format
png("nmds_plot.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

