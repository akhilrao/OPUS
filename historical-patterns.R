#####
# R script to generate figures and calculate statistics/figures of merit about historical orbital-use patterns
#####

# Load packages and source analytics-functions.R
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Read in the data UCS-JSpOC-soy-panel-22.csv, drop first column of row indices
historical_data_full = read_csv("historical-data/UCS-JSpOC-soy-panel-22.csv")[,-1]

# Filter to necessary variables
df = historical_data_full %>%
    select(COSPAR_number, Year, SATNAME_JSpOC, mean_altitude, Entity_Type, Purpose)

# Simplify Entity_Type: every row that contains Commercial should be reclassified as "Commercial", every row that does not contain "Commercial", "Starlink", or "Oneweb" should be reclassified as "Other". Then augment Entity_Type with values for Starlink and OneWeb. 
df = df %>%
    mutate(Entity_Type = ifelse(str_detect(Entity_Type, "Commercial"), "Fringe", "Other")) %>%
    mutate(Entity_Type = ifelse(str_detect(SATNAME_JSpOC, "STARLINK"), "Starlink", Entity_Type)) %>%
    mutate(Entity_Type = ifelse(str_detect(SATNAME_JSpOC, "ONEWEB"), "Oneweb", Entity_Type)) %>%
    filter(mean_altitude <= 1600)

# Construct a new variable with 35 km altitude bins starting at mean_altitude=200.
break_seq = seq(95, 1600, by = 35) # vector of breaks
label_seq = paste0(break_seq[-length(break_seq)], "-", break_seq[-1])  # vector of labels that is the same length as the break_seq
df = df %>%
    mutate(altitude_bin = cut(mean_altitude, breaks = break_seq, labels=label_seq))

# Group by year, Entity_Type, and altitude, calculate total number of satellites annually after making all Entity_Types have entries for years from 1959-2022. Then reclassify Entity_Type = Starlink or Oneweb to "Constellation". All empty or missing values of year should be filled with 0 in the satellite total. Then filter to years >= 2006
df_satshellyear = df %>%
    group_by(Year, Entity_Type, altitude_bin) %>%
    summarise(total_satellites = n()) %>%
    ungroup() %>%
    mutate(Entity_Type = ifelse(Entity_Type == "Starlink" | Entity_Type == "Oneweb", "Constellation", Entity_Type)) %>%
    complete(Year = seq(1959, 2022, 1), Entity_Type, altitude_bin, fill = list(total_satellites = 0)) %>%
    filter(Year >= 2006)

# Make heatmaps of total satellite count, year on X and altitude on Y, for each value of Entity_Type: Other, Constellation, Fringe. Use scale_fill_distiller for color. If the Entity_Type is "Fringe" use the "Blues" palette. If the Entity_Type is "Constellation" use the "Greens" palette. If the Entity_Type is "Other" use the "BuPu" palette. Then make the heatmap.

heatmap_other = df_satshellyear %>%
    filter(Entity_Type == "Other") %>%
    ggplot(aes(x = Year, y = altitude_bin, fill = total_satellites)) +
    geom_tile() +
    scale_fill_distiller(palette = "BuPu", direction = 1) +
    labs(title = "Other", x = "Year", y = "Altitude (km)", fill = "Satellites") +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme_bw()

heatmap_constellation = df_satshellyear %>%
    filter(Entity_Type == "Constellation") %>%
    ggplot(aes(x = Year, y = altitude_bin, fill = total_satellites)) +
    geom_tile() +
    scale_fill_distiller(palette = "Greens", direction = 1) +
    labs(title = "Constellations", x = "Year", y = "Altitude (km)", fill = "Satellites") +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme_bw()

heatmap_fringe = df_satshellyear %>%
    filter(Entity_Type == "Fringe") %>%
    ggplot(aes(x = Year, y = altitude_bin, fill = total_satellites)) +
    geom_tile() +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    labs(title = "Fringe", x = "Year", y = "Altitude (km)", fill = "Satellites") +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme_bw()

# Use patchwork to combine all three plots into a wide plot with collected guides at the bottom, then save each separately and the wide plot using ggsave. The joint plot plots should have width = 18, height = 6, dpi = 320, bg="white", the separate plots should have width = 8, height = 6, dpi = 320, bg="white".

ggsave("historical-data/heatmap-other.png", bg = "white", heatmap_other, width = 8, height = 10, dpi = 320)

ggsave("historical-data/heatmap-constellation.png", bg = "white", heatmap_constellation, width = 8, height = 10, dpi = 320)

ggsave("historical-data/heatmap-fringe.png", bg = "white", heatmap_fringe, width = 8, height = 10, dpi = 320)

heatmap_all = 
    (heatmap_constellation + labs(x="")) + 
    (heatmap_fringe + labs(y="") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())) + 
    (heatmap_other + labs(y="", x="") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

ggsave("historical-data/heatmap-all.png", bg = "white", heatmap_all, width = 18*0.75, height = 10*0.75, dpi = 420)

##### Compare historical fringe and constellation against 25-year rule IAM outputs

# Load analytics functions
source("analytics-functions.R")

# Adjust inputs as needed
scenarios_folder = "scenarios"
iam_25_name = "philosophical-arctic-fox-bends-LemonBalm"
iam_outputs = read_csv(paste0(scenarios_folder, "/", iam_25_name, "/", iam_25_name, "-MOCAT-equilibrium-35yrs__out-data.csv"))

# Read the files and compute altitude_bin_midpoints
altitude_bin_midpoints = read_and_compute_midpoints("x0_TLE/Counts_DEBRIS_bins_35.csv")

# Add altitude_bin_midpoints to output_data
iam_outputs_aug = data.frame(altitude_bin_midpoints, iam_outputs)

# Construct a new variable with 35 km altitude bins starting at mean_altitude=200. Doing this with a little sleight of hand in the altitude bin variable naming so that we can use the plot-heatmap function later
break_seq_2 = seq(200, 1600, by = 35) # vector of breaks
label_seq_2 = paste0(break_seq_2[-length(break_seq_2)], "-", break_seq_2[-1])  # vector of labels that is the same length as the break_seq
iam_outputs_aug = iam_outputs_aug %>%
    mutate(altitude_bin = cut(altitude_bin_midpoints, breaks = break_seq_2, labels=label_seq_2)) %>%
    select(-altitude_bin_midpoints) %>%
    rename(altitude_bin_midpoints = altitude_bin) 

# Pivot to long for plotting, filter for fringe
iam_outputs_melted = iam_outputs_aug %>%
  pivot_longer(cols = -altitude_bin_midpoints, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("type", "time"), sep = "_") %>%
  mutate(time = as.numeric(time)) %>%
  filter(type=="Su", time<=6)

# Pivot to long for plotting, filter for constellation
iam_outputs_melted_constellation = iam_outputs_aug %>%
  pivot_longer(cols = -altitude_bin_midpoints, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("type", "time"), sep = "_") %>%
  mutate(time = as.numeric(time)) %>%
  filter(type=="S", time<=6)

# Generate heatmap from model: Su
Su_plot = heatmap_plot(iam_outputs_melted, "Su", "Fringe", barwidth=1) + 
    scale_fill_distiller(name="Model", direction=1) + 
    theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()
        )

# Generate heatmap from model: S
S_plot = heatmap_plot(iam_outputs_melted_constellation, "S", "Constellation", barwidth=1) + 
    scale_fill_distiller(name="Model", direction=1, palette="YlGn") + 
    theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()
        )

# Regenerate heatmap for fringe over restricted domain
heatmap_fringe_restricted = df_satshellyear %>%
    filter(Entity_Type == "Fringe", Year >= 2017) %>%
    filter(
        altitude_bin != "95-130" & 
        altitude_bin != "130-165" &
        altitude_bin != "165-200") %>%
    ggplot(aes(x = Year, y = altitude_bin, fill = total_satellites)) +
    geom_tile() +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    labs(title = "Fringe", x = "Year", y = "Altitude (km)", fill = "Historical") +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme_bw()

# Regenerate heatmap for constellation over restricted domain
heatmap_constellation_restricted = df_satshellyear %>%
    filter(Entity_Type == "Constellation", Year >= 2017) %>%
    filter(
        altitude_bin != "95-130" & 
        altitude_bin != "130-165" &
        altitude_bin != "165-200") %>%
    ggplot(aes(x = Year, y = altitude_bin, fill = total_satellites)) +
    geom_tile() +
    scale_fill_distiller(palette = "YlGn", direction = 1) +
    labs(title = "Constellation", x = "Year", y = "Altitude (km)", fill = "Historical") +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme_bw()

# Comparison plots
comparison_plot_su = 
    (heatmap_fringe_restricted + labs(title="Fringe (historical)")) + 
    (Su_plot + labs(title="Fringe (model)") + theme_bw()) + 
    plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect')

ggsave("historical-data/comparison-plot-fringe.png", bg = "white", comparison_plot_su, width = 18*0.5*0.75, height = 10*0.75, dpi = 320)

comparison_plot_s = 
    (heatmap_constellation_restricted + labs(title="Constellation (hist)")) + 
    (S_plot + labs(title="Constellation (model)") + theme_bw()) + 
    plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect')

ggsave("historical-data/comparison-plot-constellation.png", bg = "white", comparison_plot_s, width = 18*0.5*0.75, height = 10*0.75, dpi = 320)
