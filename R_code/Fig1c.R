# this is top part of Fig 1c

# ── Load libraries ──────────────────────────────────────────────
library(readxl)    # read Excel files
library(dplyr)     # data wrangling
library(tidyr)     # pivot_longer
library(ggplot2)   # plotting

# ── 1. Read the data ────────────────────────────────────────────
df <- read_excel("Fig1b_data.xlsx")   # columns: Week, Alpha, Delta, Omicron

# ── 2. Reshape & prep labels ────────────────────────────────────
plot_df <- df %>% 
  pivot_longer(-Week, names_to = "Variant", values_to = "Count") %>% 
  mutate(
    WeekLabel = gsub("_", "_", Week),     # “wk 49 2020”
    WeekLabel = factor(WeekLabel, levels = unique(WeekLabel))  # keep order
  )

# ── 3. Define variant colours ──────────────────────────────────
variant_cols <- c(Alpha = "#B79F00",
                  Delta = "magenta",
                  Omicron = "blue")

# ── 4. Plot ────────────────────────────────────────────────────
ggplot(plot_df, aes(x = WeekLabel, y = Count, fill = Variant)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +                              # stacked bars
  scale_fill_manual(values = variant_cols, name = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of sequences") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=7),
    axis.text.y = element_text(size=8),
    axis.title.y = element_text(size=10),
    legend.position = c(0.95, 0.65),  # Adjust these coordinates to place the legend
    legend.direction = "vertical", # Make the legend horizontal
    legend.justification = c(1, 1),
    legend.text = element_text(size=10)
  )

# ── Optional: save high-resolution figure ──────────────────────
ggsave("Fig1b_plot.png", width = 4, height = 2.5, dpi = 1000)



#this is bottom part of Fig 1c
library(ggplot2)
library(readr)
library(dplyr)
library(ggthemes)
library(scales)
library(tidyverse)

# Load data
filtered_variants <- read_csv("human_trend.csv")

# Ensure date is in Date format
filtered_variants <- read_csv("human_trend.csv") %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

custom_colors <- c(
  "Alpha"   = "#B79F00",  
  "Delta"   =   "magenta",      #"#FFD900", "#FFB93B",  "#A569BD",
  "Omicron" = "#0066FF"   
)


label_data <- filtered_variants %>%
  group_by(lineage) %>%
  slice_max(order_by = proportion, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(proportion = proportion - 0.4)  # offset downward


#this chunk of code taken from Aijing's chunk of code. 
p = ggplot(filtered_variants, aes(x = date, y = proportion, colour = .data[["lineage"]],fill = .data[["lineage"]])) +
  #geom_ribbon(aes(ymin = proportion_ci_lower, ymax = proportion_ci_upper), alpha = 0.35, size = 0) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = 0, ymax = proportion), alpha = 0.4) +  # Fill the area from 0 to the proportion
  scale_x_date(date_labels = "%m/%Y", expand = c(0,0)) +
  scale_y_continuous(labels = function(x) x * 100, #scales::label_number(accuracy = 0.25)
                     expand = c(0, 0),
                     limits = c(0, 1)) +
  scale_colour_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_base()+
  labs(x = NULL, y = "Prevalence(%)") +
  geom_text(
    data = label_data,
    aes(x = date, y = proportion, label = lineage), colour = 'black',
    size = 4,  show.legend = FALSE
  ) +
  theme(legend.position = "None",
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        plot.background = element_blank(),
        plot.caption = element_text(size = 18))
p
ggsave("variant_trends_4.png", plot = p, width = 4, height = 1.5,dpi=1000)