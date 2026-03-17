# Load required libraries
library(ggplot2)
library(sf)
library(dplyr)
library(ggspatial)

# Load the saved base map
load("pa_base_map_no_X.RData")

p <- ggplot() +
  geom_sf(data = pa_counties, fill = NA, color = "black", size = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# for alpha i have C1, C2, C3 and C5 REMEMBER
# this Rscript is for C1 and C2 combined make sure to change
alpha_clust <- '(C-1, C-2)'
transmission_links_C2 <- read.csv("./C2_translinks.csv")  
transmission_links_C2$cluster <- 'C2'

transmission_links_C1 <- read.csv("./C1_translinks.csv")  
transmission_links_C1$cluster <- 'C1'

#combine the transmission links
transmission_links <- rbind(transmission_links_C1, transmission_links_C2)

# delete C2, C3... etc. too if needed
rm(transmission_links_C1)
rm(transmission_links_C2)


#location coordinates from C1 (the full one is same on all three)
location_coords <- read.delim("location_coordinates_full.tsv",
                              header = FALSE, 
                              col.names = c("location", "lat", "lon"))
# Drop deer_PA and nonPA_deer from links and coordinates
transmission_links <- transmission_links %>%
  filter(
    !FROM %in% c("deer_PA", "nonPA_deer","deer_county_unknown"),
    !TO   %in% c("deer_PA", "nonPA_deer","deer_county_unknown")
  )

location_coords <- location_coords %>%
  filter(!location %in% c("deer_PA", "nonPA_deer","deer_county_unknown"))




# Merge coordinates with transmission data
links_with_coords <- transmission_links %>%
  left_join(location_coords, by = c("FROM" = "location")) %>%
  rename(from_lat = lat, from_lon = lon) %>%
  left_join(location_coords, by = c("TO" = "location")) %>%
  rename(to_lat = lat, to_lon = lon)

links_with_coords <- links_with_coords %>%
  # keep only links with BAYES_FACTOR >= 3
  filter(BAYES_FACTOR >= 3) %>%
  mutate(
    BF_cat = case_when(
      BAYES_FACTOR >= 3    & BAYES_FACTOR <= 10   ~ "3–10",
      BAYES_FACTOR >  10   & BAYES_FACTOR <= 100  ~ "10–100",
      BAYES_FACTOR >  100  & BAYES_FACTOR <= 1000 ~ "100–1000",
      BAYES_FACTOR >  1000                        ~ ">1000"
    ),
    BF_cat = factor(BF_cat,
                    levels = c("3–10", "10–100", "100–1000", ">1000")),
    cluster_label = paste0("C-", sub("^C", "", cluster))
  )
bf_levels <- c("3–10", "10–100")

dummy_links <- data.frame(
  from_lon = NA_real_,
  from_lat = NA_real_,
  to_lon   = NA_real_,
  to_lat   = NA_real_,
  BF_cat   = factor(bf_levels, levels = bf_levels)
)
# Create the transmission map by adding arrows to the base plot
transmission_map <- p +
  # Move elevation legend downward
  theme(
    legend.position = c(0.55, 0.95),  # This positions the main legend box
    legend.direction = "horizontal",  # Arrange legends horizontally
    legend.box = "horizontal",         # Box arrangement
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title.position = 'top'
  ) +
  
  # Add transmission links as curved arrows with improved styling
  geom_curve(
    data = links_with_coords,
    aes(x = from_lon, y = from_lat, 
        xend = to_lon, yend = to_lat,
        linewidth = BF_cat,
        colour = cluster_label),
    curvature = 0.28,
    arrow = arrow(
      length = unit(0.25, "cm"),
      type = "open",
      angle = 30,
      ends = "last"
    ),
    alpha = 0.8,
    lineend = "round",
    show.legend = c(colour = TRUE, linewidth = FALSE)
  ) +
  geom_curve(
    data = dummy_links,
    aes(x = from_lon, y = from_lat, 
        xend = to_lon, yend = to_lat, linewidth = BF_cat),
    curvature = 0.28,
    arrow = arrow(
      length = unit(0.25, "cm"),
      type = "open",
      angle = 30,
      ends = "last"
    ),
    lineend = "round",
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +

  annotate(
    "text",
    x = -80,           # adjust
    y = 41.95,          # adjust
    label = 'Alpha',   # your text
    size = 4
  ) +
 # annotate(
#    "text",
#    x = -79.88,           # adjust
 #   y = 41.75,          # adjust
  #  label = alpha_clust,
   # size = 3.5
  #) +
  scale_linewidth_manual(
    name   = "Bayes factor",
    values = c(
      "3–10"     = 0.4,
      "10–100"   = 0.9),
    breaks = c("3–10", "10–100"),
    drop   = FALSE
  ) +
  scale_color_manual(
    name   = "Cluster",
    values = c("C-1" = "#B79F00", "C-2" = "black"),
    breaks = c("C-1", "C-2")
  ) +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      keywidth = 0.6, keyheight = 0.3,
      default.unit = "lines",
      order = 1,
      override.aes = list(
        colour = c("#B79F00", "black"),
        alpha = 1,
        linewidth = 0.4
      )
    ),
    linewidth = guide_legend(title.position = 'top', title.hjust = 0.5,
          keywidth = 0.6, keyheight = 0.3, default.unit = 'lines'), fill = 'none'
  )

# Display the map
print(transmission_map)

# Save the final map
ggsave(paste0("phylogeography4_no_outside_PA_v2_",alpha_clust,".png"), plot = transmission_map, 
       width = 3.7, height = 2.7, dpi = 700)