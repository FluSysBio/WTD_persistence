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

# Read your transmission data

#this one is D5 and D6 combined, doubcle check transmission link variable
delta_clust <- '(D-5, D-6)'

transmission_links_D5 <- read.csv("D5_translinks.csv")%>% mutate(cluster = "D-5")  
transmission_links_D6 <- read.csv("D6_translinks.csv")%>% mutate(cluster = "D-6")  

transmission_links <- rbind(transmission_links_D5,transmission_links_D6)

rm(transmission_links_D5)
rm(transmission_links_D6)

location_coords <- read.delim("locationCoordinates_full.txt", header = FALSE, 
                              col.names = c("county", "lat", "lon"))

# Drop deer_PA and nonPA_deer from links and coordinates
transmission_links <- transmission_links %>%
  filter(
    !FROM %in% c("deer_PA", "nonPA_deer","deer_county_unknown","Human"),
    !TO   %in% c("deer_PA", "nonPA_deer","deer_county_unknown","Human")
  )

location_coords <- location_coords %>%
  filter(!county %in% c("deer_PA", "nonPA_deer","deer_county_unknown","Human"))

# adjusting coordinates of Delaware county for better arrow positioning
location_coords <- location_coords %>%
  mutate(lat = if_else(county == "Delaware", 39.875473, lat),
         lon = if_else(county == "Delaware", -75.421619, lon))


# Merge coordinates with transmission data
links_with_coords <- transmission_links %>%
  left_join(location_coords, by = c("FROM" = "county")) %>%
  rename(from_lat = lat, from_lon = lon) %>%
  left_join(location_coords, by = c("TO" = "county")) %>%
  rename(to_lat = lat, to_lon = lon)

links_with_coords <- links_with_coords %>%
  mutate(
    from_lon = as.numeric(from_lon),
    to_lon   = as.numeric(to_lon),
    from_lat = as.numeric(from_lat),
    to_lat   = as.numeric(to_lat)
  )

links_with_coords <- links_with_coords %>%
  # keep only links with BAYES_FACTOR >= 3, same as code block A
  filter(BAYES_FACTOR >= 3) %>%
  mutate(
    BF_cat = case_when(
      BAYES_FACTOR >= 3    & BAYES_FACTOR <= 10   ~ "3–10",
      BAYES_FACTOR >  10   & BAYES_FACTOR <= 100  ~ "10–100",
      BAYES_FACTOR >  100  & BAYES_FACTOR <= 1000 ~ "100–1000",
      BAYES_FACTOR >  1000                        ~ ">1000"
    ),
    BF_cat = factor(BF_cat,
                    levels = c("3–10", "10–100", "100–1000", ">1000"))
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
        color = cluster),
    curvature = 0.28,
    arrow = arrow(
      length = unit(0.25, "cm"),
      type = "open",
      angle = 30,
      ends = "last"
    ),
    alpha = 0.8,
    lineend = "round"
  ) +
  annotate(
    "text",
    x = -80,           # adjust
    y = 41.95,          # adjust
    label = "Delta",   # your text
    size = 4
  ) +
  #annotate(
  #  "text",
  #  x = -79.88,           # adjust
  #  y = 41.75,          # adjust
  #  label = delta_clust,
  #  size = 3.5
  #) +
  scale_linewidth_manual(
    values = c(
      "3–10"     = 0.4,
      "10–100"   = 0.9
    ),
    breaks = c("3–10", "10–100"),
    drop   = FALSE
  ) +
  scale_color_manual(
    values = c("D-5" = "black", "D-6" = "magenta"),
    breaks = c("D-5", "D-6"),
    name   = "Cluster"
  ) +
  guides(
    linewidth = 'none',
    fill = 'none',
    color = guide_legend(title.position = "top",
                         title.hjust = 0.5,
                         keywidth = 0.6, keyheight = 0.3,
                         default.unit = "lines",
                         override.aes = list(linewidth = 0.4, alpha = 1))
  )

# Display the map
print(transmission_map)

# Save the final map
ggsave(paste0("phylogeography_v7_",delta_clust,"_v2.png"), plot = transmission_map, 
       width = 3.7, height = 2.7, dpi = 700)