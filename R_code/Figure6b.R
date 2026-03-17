
library(ggplot2)
library(reshape2)
library(geosphere)
library(dplyr)
library(ecodist)
library(readxl)

patristic_matrix <- as.matrix(read.csv("pairwise_patristic_distances_delta.csv",
                                       row.names = 1, check.names = FALSE))

# Fix the naming issue: replace spaces with underscores in the patristic matrix names
rownames(patristic_matrix) <- gsub(" ", "_", rownames(patristic_matrix))
colnames(patristic_matrix) <- gsub(" ", "_", colnames(patristic_matrix))

# Read the sequence-county mapping file
sequence_county <- read.csv("./sequence_list_county_colors_updated.csv",
                            stringsAsFactors = FALSE)


# Read the county coordinates data
county_coords <- read.delim("./location_coordinates.tsv", 
                            stringsAsFactors = FALSE,header = FALSE, 
                            col.names = c('county','Latitude','Longitude'))

# Filter the patristic matrix to include only samples with valid counties
valid_samples <- sequence_county$Name

common_samples <- intersect(rownames(patristic_matrix), valid_samples)

# Check if we now have matches
if(length(common_samples) == 0) {
  stop("Still no common samples found after fixing naming issues. Please check manually.")
}

print(paste("Found", length(common_samples), "common samples"))

# Read cluster labels (Excel) and clean labels the same way (spaces -> underscores)
cluster_df <- readxl::read_xlsx("./delta_PA_deer_clusterlabels.xlsx")

# Keep only D6 samples
D6_labels <- cluster_df$label[cluster_df$Cluster == "D-6"]

# Intersect with common_samples to ensure they exist in both patristic matrix and county list
D6_samples <- intersect(common_samples, D6_labels)

if (length(D6_samples) == 0) {
  stop("No overlap between D6 cluster labels and common_samples.")
}

# From here on, only use D6 samples
common_samples <- D6_samples
print(paste("Using", length(common_samples), "D6 samples"))

patristic_filtered <- patristic_matrix[common_samples, common_samples]

# Create a mapping of sample to county coordinates
sample_coords <- inner_join(sequence_county, county_coords, 
                            by = c("County" = "county"))
sample_coords <- sample_coords[sample_coords$Name %in% common_samples, ]
rownames(sample_coords) <- sample_coords$Name

# Ensure the sample order matches the patristic matrix order
sample_coords <- sample_coords[common_samples, ]

sample_coords <- sample_coords %>%
  mutate(
    Longitude = as.numeric(as.character(Longitude)),
    Latitude = as.numeric(as.character(Latitude))
  )

# Calculate geographic distance matrix
geo_dist_matrix <- distm(sample_coords[, c("Longitude", "Latitude")],
                         fun = distHaversine) / 1000
rownames(geo_dist_matrix) <- common_samples
colnames(geo_dist_matrix) <- common_samples

# Convert to dist objects for Mantel test
patristic_dist <- as.dist(patristic_filtered)
geo_dist <- as.dist(geo_dist_matrix)

# Perform Mantel test with bootstrapping using ecodist
set.seed(123)
mantel_result_ecodist <- ecodist::mantel(patristic_dist ~ geo_dist, 
                                         nperm = 9999, nboot = 5000,
                                         mrank = FALSE)

# The ecodist mantel function doesn't directly provide CI, but we can calculate from permutations
# Extract the Mantel statistic and permutation distribution
mantel_r <- mantel_result_ecodist['mantelr']
permutation_values <- mantel_result_ecodist[2:length(mantel_result_ecodist)]

# Calculate confidence intervals
#ci_95 <- quantile(permutation_values, probs = c(0.025, 0.975))
#ci_90 <- quantile(permutation_values, probs = c(0.05, 0.95))
#ci_75 <- quantile(permutation_values, probs = c(0.125, 0.875))
lower_ci <- mantel_result_ecodist['llim.2.5%']
upper_ci <- mantel_result_ecodist['ulim.97.5%']

# Print results
print(paste("Mantel r =", round(mantel_r, 3)))
#print(paste("95% CI: [", round(ci_95[1], 3), ",", round(ci_95[2], 3), "]"))
#print(paste("90% CI: [", round(ci_90[1], 3), ",", round(ci_90[2], 3), "]"))
#print(paste("75% CI: [", round(ci_75[1], 3), ",", round(ci_75[2], 3), "]"))

# Calculate p-value (proportion of permutations with r >= observed r)
#p_value <- sum(permutation_values >= mantel_r) / length(permutation_values)
p_value <- mantel_result_ecodist['pval3']
print(paste("p-value =", p_value))
cat(sprintf("95%% Confidence Interval: (%.3f, %.3f)\n", lower_ci, upper_ci))
# Prepare data for plotting
plot_data <- data.frame(
  Geographic = as.vector(geo_dist_matrix[upper.tri(geo_dist_matrix)]),
  Genetic = as.vector(patristic_filtered[upper.tri(patristic_filtered)])
)

# Create scatter plot with confidence interval in subtitle
gen_geo_plot <- ggplot(plot_data, aes(x = Geographic, y = Genetic)) +
  geom_point(alpha = 0.6, size = 0.5, color = 'magenta', shape = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", 
              linewidth = 0.5, fill = "grey80") +
  labs(
    x = "Geographic Distance (km)",
    y = "Patristic Distance",
    #title = "Delta: D-6",
    subtitle = sprintf("Mantel r = %.3f, p = %.4f",
                       mantel_r, p_value)
  ) +
  coord_cartesian(xlim = c(0, 250),   # <-- set your limits here
                  ylim = c(0, 0.007)) +
  theme_classic() +
  theme(#plot.title = element_text(hjust = 0.85, vjust=-8,face = "bold",size = 9),
        axis.title.x = element_text(size = 8),
        #axis.title.y = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 7),
        #axis.text.y = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.85, vjust=-6, size=8)) 
  #annotate(
  #  "text",
  #  x = 200,
  #  y = 0.001,
  #  label = 'Note: x-axis limit',
  #  size = 3
  #)

# Display and save the plot
print(gen_geo_plot)
write.csv(plot_data,'Mantel_test_geographic_delta_D6.csv' ,row.names=FALSE)
ggsave("Mantel_test_geographic_delta_D6.png",gen_geo_plot,
       width = 3, height = 2, dpi = 600)
