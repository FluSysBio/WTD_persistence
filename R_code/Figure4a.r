#this is Figure 4a

library(scales)
library(ggplot2)
library(ggtree)
library(phytools)
library(treeio)
library(TreeTools)
library(ape)
library(tidytree)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

options(ignore.negative.edge=TRUE)


tree <- read.beast('alpha_tree2_dedup.ccv_output_selected.trees.output_58')

#TEST = read.csv('alpha_deer_date.txt',sep = "\t",header = FALSE)
#TEST = read.csv("../../../transmission_rate/WTD_clusters.csv",sep = ",")
#TEST = read.csv("tip_names_in_plot_order_addCluster.csv",sep = ",")
TEST = read.csv('TableS5_alphaCluster_OLD_v2_single_human_and_C9_corrected.csv')

#NOTE: keep in mind the color scheme and cluster labeling was manual. 
#The figures powerpoint file as of Oct 2025 has the up to date cluster label scheme

default_colors <- scales::hue_pal()(12)  

color_mapping <- c(
  "Human" = "lightgrey",
  "Other WTD" = "black",
  '1' = default_colors[1],
  '2' = default_colors[2],
  '3' = default_colors[3],
  '4' = default_colors[4],
  '5' = default_colors[5],
  '6' = default_colors[6],
  '7' = default_colors[7],
  '8' = default_colors[8],
  '9' = default_colors[9],
  '10' = default_colors[10],
  '11' = default_colors[11],
  '12' = default_colors[12],
  'X' = default_colors[11]
)


TEST$cluster=factor(TEST$cluster, 
                    levels = names(color_mapping),
                    ordered = TRUE)
TEST$Host2 <- ifelse(TEST$cluster == "Human", "Human", "WTD")
p <- ggtree(tree,mrsd="2023-03-15",size=0.5) %<+% TEST +geom_treescale(x=2022.3,y=5,offset = 3,color='black',width=0.5)+
  #geom_tiplab(aes(color=factor(cluster)))+
  #geom_tippoint(aes(color=factor(cluster)), shape = 1, size = 1.5,
  #              alpha = 0.95,
  #              stroke = 0.3)+
  theme_tree2()+
  theme(legend.position = c(0.20,0.97), 
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 8,colour = 'black'),
        plot.margin   = margin(t = 5.5, r = 20, b = 10, l = 20, unit = "pt")) +
  #geom_text2(color = "Magenta",aes(subset=!isTip, label=node), hjust=-0.3,size=3)+
  #geom_nodelab(color = "Orange",aes(x=branch, label=round(as.double(posterior),2), subset=posterior>0.69999), vjust=0,hjust=0.5, size=3)+
  scale_color_manual(
    values = c("Human" = "lightgrey",
               "WTD"  = default_colors[3]),
    breaks = c("Human", "WTD"),
    name = "",
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        fill  = c("lightgrey", default_colors[3])  
      )
    )
  )+
  scale_fill_manual(
    values = c("Human" = "lightgrey",
               "WTD"   = default_colors[3]),
    guide = "none"
  ) +
  annotate(
    "text",
    x     = 2019.5,     # time axis in same units as tree
    y     = 310,         # tip index / vertical position
    label = "Lineage",
    hjust = 0,
    size  = 3.528
  )+
  annotate(
    "text",
    x     = 2019.7,     # time axis in same units as tree
    y     = 295,         # tip index / vertical position
    label = "B*",
    hjust = 0,
    size  = 2.864
  )

#shifting tips a bit to right
plot_data <- p$data

x_shift_tips <- max(plot_data$x, na.rm = TRUE) * 0.00002
shifted_tips <- plot_data %>%
  filter(isTip) %>%
  mutate(x = x + x_shift_tips)
p <- p + geom_tippoint(
  data = shifted_tips,
  aes(
    x = x,
    y = y,
    color = Host2,
    fill = Host2
  ), shape = 21, size = 1.5,
  alpha = 0.95,
  stroke = 0.3,
  inherit.aes = FALSE
)

# add vertical padding so the top tip isn't cut off
y_min <- min(plot_data$y, na.rm = TRUE)
y_max <- max(plot_data$y, na.rm = TRUE)
p <- p + coord_cartesian(ylim = c(y_min, y_max + 5))

print(p)


tip_csv <- plot_data %>%
  dplyr::filter(isTip) %>%
  dplyr::arrange(dplyr::desc(y)) %>%   # top-to-bottom; use arrange(y) for bottom-to-top
  dplyr::transmute(label) 

#write.csv(tip_csv,
#          file = "tip_labels_in_plot_order_alpha_tree_Aijing.csv",
#          row.names = FALSE,
#          quote = TRUE)

#ggsave("./clade__tree.pdf", width = 200, height = 200, units = "cm", limitsize = FALSE)
ggsave("alpha_tree_v2.png", p, width = 3.5, height = 6, dpi = 600,limitsize=FALSE)

p<-ggtree(tree,mrsd="2023-03-15",size=0.7)
tip_names_in_plot_order<-get_taxa_name(p)

print(p)
#write.csv(tip_names_in_plot_order, 'tip_names_in_plot_order.csv', row.names = FALSE)

