
library(dplyr)
library(ggplot2)
library(reshape2)
library(geosphere)
library(ecodist)

patristic_matrix <- as.matrix(read.csv("pairwise_patristic_distances.csv",
                                       row.names = 1, check.names = FALSE))

# Fix the naming issue: replace spaces with underscores in the patristic matrix names
rownames(patristic_matrix) <- gsub(" ", "_", rownames(patristic_matrix))
colnames(patristic_matrix) <- gsub(" ", "_", colnames(patristic_matrix))

# Read the sequence-county mapping file
sequence_county <- read.delim(".sequence_county.tsv", stringsAsFactors = FALSE)

cluster_info <- read.delim("./sequence_cluster_CORRECT.tsv", stringsAsFactors = FALSE)
cluster_info <- as_tibble(cluster_info)
cluster_info <- cluster_info[!cluster_info$cluster %in% c('HomoSapiens','Other WTD'),]


# Filter out samples with imprecise locations (human, nonPA_deer, deer_PA)
valid_counties <- sequence_county[!sequence_county$county %in% c("human", "nonPA_deer", "deer_PA"), ]

# Read the county coordinates data
county_coords <- read.delim("./combined_rates/phylo1/location_coordinates_No_human.tsv", 
                            stringsAsFactors = FALSE,header = FALSE, 
                            col.names = c('county','Latitude','Longitude'))

# Filter the patristic matrix to include only samples with valid counties
valid_samples <- valid_counties$Name

common_samples <- intersect(rownames(patristic_matrix), valid_samples)

# Check if we now have matches
if(length(common_samples) == 0) {
  stop("Still no common samples found after fixing naming issues. Please check manually.")
}

print(paste("Found", length(common_samples), "common samples"))

patristic_filtered <- patristic_matrix[common_samples, common_samples]

# Create a mapping of sample to county coordinates
sample_coords <- merge(valid_counties, county_coords, by.x = "county", by.y = "county")
sample_coords <- sample_coords[sample_coords$Name %in% common_samples, ]
rownames(sample_coords) <- sample_coords$Name

# Ensure the sample order matches the patristic matrix order
sample_coords <- sample_coords[common_samples, ]

# Calculate geographic distance matrix
geo_dist_matrix <- distm(sample_coords[, c("Longitude", "Latitude")], fun = distHaversine) / 1000
rownames(geo_dist_matrix) <- common_samples
colnames(geo_dist_matrix) <- common_samples

message("Make sure check out this cluster number here.")
# ===== Subset to a chosen cluster =====
selected_cluster <- 1   

# Samples that belong to the chosen cluster AND are in common_samples
cluster_samples <- cluster_info$seqName[cluster_info$cluster == selected_cluster]
cluster_samples <- intersect(cluster_samples, common_samples)

if (length(cluster_samples) < 2) {
  stop(sprintf("Not enough samples in cluster %s after intersection (need >= 2).", selected_cluster))
}

# Subset the square matrices to these samples
patristic_sub   <- patristic_filtered[cluster_samples, cluster_samples]
geo_dist_sub    <- geo_dist_matrix[cluster_samples, cluster_samples]

# Convert to dist objects for Mantel test
patristic_dist <- as.dist(patristic_sub)
geo_dist <- as.dist(geo_dist_sub)

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
  Geographic = as.vector(geo_dist_sub[upper.tri(geo_dist_sub)]),
  Genetic    = as.vector(patristic_sub[upper.tri(patristic_sub)])
)



# Create scatter plot with confidence interval in subtitle
gen_geo_plot <- ggplot(plot_data, aes(x = Geographic, y = Genetic)) +
  geom_point(color="#B79F00",alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", 
              linewidth = 0.5, fill = "grey80") +
  labs(
    x = "Geographic Distance (km)",
    y = "Patristic Distance",
    #title = sprintf("Alpha: C-%s", selected_cluster),
    subtitle = sprintf("Mantel r = %.3f, p = %.4f",
                       mantel_r,  p_value)
  ) +
  coord_cartesian(xlim = c(0, 250),   # <-- set your limits here
                  ylim = c(0, 0.007)) +
  theme_classic() +
  theme(#plot.title = element_text(hjust = 0.85, vjust=-8,face = "bold",size = 9),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.85, vjust=-6, size=8))

# Display and save the plot
print(gen_geo_plot)
write.csv(plot_data, sprintf("geographic_vs_genetic_Mantel_C-%s.csv",selected_cluster),
          row.names = FALSE)

ggsave(sprintf("geographic_vs_genetic_Mantel_C-%s.png",selected_cluster),
       plot = gen_geo_plot, width = 3, height = 2, dpi = 600)