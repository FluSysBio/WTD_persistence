# Load necessary libraries
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(lubridate)
library(readxl)

mcc_tree <- read.beast("./combined_tree1_treeAnnotator")
date_file <- "../sequence_dates.tsv"
host_file <- '../sequence_host.tsv'


host_info <- read.table(host_file, header = TRUE, sep = '\t')
cluster_df <- read_xlsx("delta_PA_deer_clusterlabels.xlsx")
lineage_file <- "extracted_304_deer_human_sequences_from_05_PANGO.csv"
lineage_df <- read.csv(lineage_file,
                       stringsAsFactors = FALSE)
tip_dates    <- read.table(date_file,  header = TRUE, sep = "\t")
tip_dates <- tip_dates %>%
  mutate(decimal_date = decimal_date(ymd(date)))

tip_dates <- tip_dates %>%
  mutate(date = ymd(date))  # keep decimal_date if you still need it elsewhere

mrsd_date <- max(tip_dates$date, na.rm = TRUE)
message('All Delta deer are assumed to be from Pennsylvania')
# Attach Host and date info to the tree tips
tree_df <- mcc_tree %>%
  as_tibble() %>%
  left_join(host_info, by = c("label" = "Name")) %>%    # adds 'Host'
  left_join(tip_dates,  by = c("label" = "header")) %>% # adds 'date' 
  left_join(cluster_df, by = 'label') %>%
  left_join(lineage_df, by = c('label' = 'Sequence.name'))

# Plot: circles at tips colored by Host, time scale at bottom
p <- ggtree(mcc_tree_with_data,
            mrsd = mrsd_date,
            ladderize = TRUE,
            right = TRUE) +
  geom_tippoint(
    aes(
      color = Host
    ),
    size = 1.5,
    shape = 16,
    alpha = 0.95,
    stroke = 0.3,
    show.legend = TRUE
  ) +
  scale_color_manual(
    # This is STRICTLY kept as requested (controls the point OUTLINE)
    values = c("Human" = "lightgrey", "WTD" = "magenta"),
    labels = c('Human' = "Human",
               "WTD"   = "WTD"),
    na.value = "blue",
    name = ""
  ) +
  theme_tree2() +                       # adds the time axis at the bottom
  scale_x_continuous(
    limits = c(2019, 2023),
    breaks = c(2019, 2020, 2021, 2022, 2023),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  scale_y_reverse() +
  theme(
    plot.title    = element_text(hjust = 0.5),
    legend.position = c(0.35, 0.95),
    legend.title  = element_text(size = 10),
    legend.text   = element_text(size = 10),
    legend.box    = "vertical",
    axis.text.x   = element_text(size = 8, color='black'),
    plot.margin   = margin(t = 5.5, r = 20, b = 10, l = 20, unit = "pt")
  )+
  guides(
    color = guide_legend(order = 1)
  )
# Get plotted data from ggtree
plot_data <- p$data

# Pick one WTD tip per cluster (excluding singletons "S")
x_shift <- max(plot_data$x, na.rm = TRUE) * 0.00004

# D-clusters: one representative per D1, D2, ...
cluster_ann_D <- plot_data %>%
  dplyr::filter(
    isTip,
    Host == "WTD",
    !is.na(Cluster),
    Cluster != "S"         # exclude singletons here
  ) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::slice(1) %>%         # one tip per D-cluster
  dplyr::ungroup() %>%
  dplyr::mutate(x = x + x_shift)


p <- p +
  geom_text(
    data  = cluster_ann_D,
    aes(x = x, y = y, label = Cluster),
    color = "black",
    hjust = 0,
    vjust = 0.5,
    size  = 3
  ) 
print(p)
ggsave("delta_tree_clusterlabelwith_lineageshape_v3.png",
       p, width = 3.5, height = 6, dpi = 600,limitsize=FALSE)