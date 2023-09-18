#####
# R script to generate figures and calculate statistics/figures of merit about the simulation. 
# "Constellation" and "Fringe" are used as synonyms for "Constellation" and "Fringe" respectively.
#####

# Load packages and source analytics-functions.R
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(data.table)
source("analytics-functions.R")

# Read in command line arguments: stem, model_type, and launch_pattern_type, and model_horizon
# args = commandArgs(trailingOnly = TRUE)
# stem = args[1]
# model_type = args[2]
# launch_pattern_type = args[3]
# model_horizon = args[4]

# Script version
stem = "staying-coral-stew-Soybean"
model_type = "GMPHD"
launch_pattern_type = "equilibrium"
model_horizon = "10"

# Construct input filename
input_name = paste0(stem, "-", model_type, "-", launch_pattern_type, "-", model_horizon, "yrs") #Input filename

# Read the files and compute altitude_bin_midpoints
altitude_bin_midpoints = read_and_compute_midpoints("x0_TLE/Counts_DEBRIS_bins_35.csv")

# Read the output CSV
output_data = read_csv(paste0("scenarios/", stem, "/", input_name, "__out-data.csv"))

# # Read the collision parameters CSV
# collision_parameters = cbind(
#   altitude_bin_midpoints,
#   read_csv(paste0("scenarios/", stem, "/", stem, "--collision-parameters.csv"))
#   )

# # Read the fragmentation parameters CSV
# fragmentation_parameters = read_csv(paste0("scenarios/", stem, "/fragmentation-parameters.csv"))

# Add altitude_bin_midpoints to output_data
output_data_aug = data.frame(altitude_bin_midpoints, output_data)

# Reshape data for ggplot
output_data_melted = output_data_aug %>%
  pivot_longer(cols = -altitude_bin_midpoints, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("type", "time"), sep = "_") %>%
  mutate(time = as.numeric(time))

# Attach expected loss as new column
output_data_melted = output_data_melted %>% filter(type=="S" | type=="Su" | type=="Dlt" | type=="Dst" | type=="Dlnt")

# Compute total collision_probability for each type and each across all altitude_bin_midpoints
aggregate_metrics_paths = output_data_melted %>%
  filter(type=="S" | type=="Su" | type=="Dlt" | type=="Dst" | type=="Dlnt") %>% 
  group_by(type, time) %>%
  summarise(
    total_objects = sum(value, na.rm=TRUE)
    )
output_data_melted %>%
  filter(type=="Dst")

aggregate_metrics_paths %>% print(n=Inf)

write_csv(aggregate_metrics_paths, paste0("scenarios/", stem, "/", input_name, "__aggregate-metrics-paths.csv"))

#####
# Generate plots
#####

# Plot total collision probability for each type over time. Replace "type" with "Type" in the legend, and replace labels as D=Derelict, S=Constellation, Su=Fringe.
## Extract the first three colors from the Set2 palette
color_set <- brewer.pal(3, "Set1")
## Generate plot

# Make and combine heatmap plots
plot_Su = output_data_melted %>%
  filter(type == "Su") %>%
  ggplot(aes(x = time, y = altitude_bin_midpoints, fill = value)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Blues", 
    guide = guide_colorbar(barwidth = 8), 
    direction=1,
    name="Fringe satellites"
    ) +
  labs(y="Altitude[km]", x="Time", title="Fringe satellites") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme_bw()

plot_Dlnt = output_data_melted %>%
  filter(type == "Dlnt") %>%
  ggplot(aes(x = time, y = altitude_bin_midpoints, fill = value)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "OrRd", 
    guide = guide_colorbar(barwidth = 8), 
    direction=1,
    name="LNT"
    ) +
  labs(y="Altitude[km]", x="Time", title="Lethal non-trackable debris") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme_bw()

plot_Dst = output_data_melted %>%
  filter(type == "Dst") %>%
  ggplot(aes(x = time, y = altitude_bin_midpoints, fill = value)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "OrRd", 
    guide = guide_colorbar(barwidth = 8), 
    direction=1,
    name="Small trackables"
    ) +
  labs(y="Altitude[km]", x="Time", title="Small trackable debris") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme_bw()

plot_Dlt = output_data_melted %>%
  filter(type == "Dlt") %>%
  ggplot(aes(x = time, y = altitude_bin_midpoints, fill = value)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "OrRd", 
    guide = guide_colorbar(barwidth = 8), 
    direction=1,
    name="Large trackables"
    ) +
  labs(y="Altitude[km]", x="Time", title="Large trackable debris") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  theme_bw()

final_plot = ((plot_Dlt / plot_Dst) | (plot_Su / plot_Dlnt)) + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
final_plot = final_plot & theme(legend.position = "bottom")

ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots.png"), final_plot, width = 12, height = 8, dpi = 320, bg="white")
