#####
# R script to generate figures and calculate statistics/figures of merit about the simulation. 
# "Slotted" and "Unslotted" are used as synonyms for "Constellation" and "Fringe" respectively.
#####

### MAIN SCRIPT BLOCK ###

# Load packages
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

source("analytics-functions.R")

# Read in command line arguments: stem, model_type, and launch_pattern_type, and model_horizon
args = commandArgs(trailingOnly = TRUE)
stem = args[1]
model_type = args[2]
launch_pattern_type = args[3]
model_horizon = args[4]
out_folder = args[5]

stem_1 = args[1]
stem_2 = args[2] 

model_type_1 = args[3]
model_type_2 = args[4]

launch_pattern_type_1 = args[5]
launch_pattern_type_2 = args[6]

model_horizon_1 = args[7]
model_horizon_2 = args[8]

scenario_label_1 = args[9]
scenario_label_2 = args[10]

out_folder = args[11]

# Construct input filename
# input_name = paste0(stem, "-", model_type, "-", launch_pattern_type, "-", model_horizon, "yrs") #Input filename

# Script version
# stem_1 = "polluting-flamingo-enunciate-Cilantro"
# model_type_1 = "MOCAT"
# launch_pattern_type_1 = "sat-feedback"
# model_horizon_1 = "35"
# input_name_1 = "polluting-flamingo-enunciate-Cilantro-MOCAT-sat_feedback-35yrs"
# scenario_label_1 = "sat-feedback"

# stem_2 = "polluting-flamingo-enunciate-Cilantro"
# model_type_2 = "MOCAT"
# launch_pattern_type_2 = "equilibrium"
# model_horizon_2 = "35"
# input_name_2 = "polluting-flamingo-enunciate-Cilantro-MOCAT-equilibrium-35yrs"
# scenario_label_2 = "equilibrium"

# out_folder = stem_1
# # out_folder = "5 vs 25 year eol comparison"
# # out_folder = "5 year no-ouf vs 25 year ouf comparison"
        
# stem_1 = "polluting-flamingo-enunciate-Cilantro" 
# stem_2 = "polluting-flamingo-enunciate-Cilantro" 

# model_type_1 = "MOCAT"
# model_type_2 = "MOCAT"

# launch_pattern_type_1 = "sat_feedback"
# launch_pattern_type_2 = "equilibrium"

# model_horizon_1 = "35"
# model_horizon_2 = "35"

# scenario_label_1 = "Satellite feedback behavior"
# scenario_label_2 = "Open-access behavior" 

# out_folder = "comparison--benchmark-eqm-satfeedback"

input_name_1 = paste0(stem_1, "-", model_type_1, "-", launch_pattern_type_1, "-", model_horizon_1, "yrs") #Input filename
input_name_2 = paste0(stem_2, "-", model_type_2, "-", launch_pattern_type_2, "-", model_horizon_2, "yrs") #Input filename

print(input_name_1)
print(input_name_2)

# Read the files and compute altitude_bin_midpoints
altitude_bin_midpoints = read_and_compute_midpoints("x0_TLE/Counts_DEBRIS_bins_35.csv")

## Use if you want to use metrics already computed in analytics.R. Faster, more consistent.
comparison_data_1 = read_csv(paste0("scenarios/", stem_1, "/", input_name_1, "__aggregate-metrics-paths.csv")) 
# %>%
#   group_by(type, time) 
  # %>%
  # mutate(normalized_ssr_index = ssr_index/ssr_index[1])
comparison_data_2 = read_csv(paste0("scenarios/", stem_2, "/", input_name_2, "__aggregate-metrics-paths.csv"))

# Join the two dataframes
aggregate_metrics_comparison_wide = left_join(comparison_data_1, comparison_data_2, by=c("type", "time"), suffix=c("__1", "__2"))
# Pivot longer
aggregate_metrics_comparison_long = pivot_longer(
  aggregate_metrics_comparison_wide, 
  cols = -c("type", "time"), 
  names_to = "variable", 
  values_to = "values") %>%
  separate(variable, into = c("variable", "label"), sep = "__") %>%
  mutate(
    label = ifelse(label == "1", scenario_label_1, scenario_label_2),
    type = case_when(
      type == "D" ~ "Derelict",
      type == "Su" ~ "Fringe",
      type == "S" ~ "Constellation"
    )
    )

aggregate_metrics_comparison_long
aggregate_metrics_comparison_long %>% filter(type=="Fringe")

#####
# Generate plots
#####

# Plot ecob index for each type over time. Replace "type" with "Type" in the legend, and replace labels as D=Derelict, S=Slotted, Su=Unslotted.
## Function to generate comparison plots. Currently only does ssr index, but can do others.
comp_plot_fn <- function(df, label_type="Satellite feedback behavior", time_window=20, plot_title="Open-access", variable_name="ssr_index", leftmost=FALSE){

# Extract the first three colors from the Set2 palette
color_set <- brewer.pal(3, "Set1")

# df=aggregate_metrics_comparison_long
# label_type="Satellite feedback behavior"
# time_window=20
# plot_title="Open-access"
# variable_name="ssr_index"
# leftmost=FALSE

  if(variable_name=="expected_max_loss") {
    plot_base = df %>%
      filter(variable==variable_name, label==label_type, type=="Fringe") %>%
    ggplot(aes(x = time)) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Expected maximum welfare [M$]"
    ) +
    scale_y_continuous(labels = label_number(scale = 1e-6, big.mark = ","))
  }

  if(variable_name=="ssr_index") {
    plot_base = df %>%
      filter(variable==variable_name, label==label_type, time<=time_window) %>%
    ggplot(aes(x = time, color = type)) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Aggregate SSR index"
    ) +
    scale_color_manual(
      name = "Object type",
      breaks = c("Derelict", "Constellation", "Fringe"),
      labels = c("Derelict", "Constellation", "Fringe"),
      values = c("Derelict" = color_set[1], "Constellation" = color_set[2], "Fringe" = color_set[3])
    )
  }

  if(variable_name=="normalized_ssr_index") {
    plot_base = df %>%
      filter(variable==variable_name, label==label_type, time<=time_window) %>%
    ggplot(aes(x = time, color = type)) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Normalized aggregate SSR index",
      color = "Object type"
    ) +
    scale_color_manual(
      name = "Object type",
      breaks = c("Derelict", "Constellation", "Fringe"),
      labels = c("Derelict", "Constellation", "Fringe"),
      values = c("Derelict" = color_set[1], "Constellation" = color_set[2], "Fringe" = color_set[3])
    )
  }

  if(variable_name=="total_objects") {
    plot_base = df %>%
      filter(variable==variable_name, label==label_type, time<=time_window) %>%
    ggplot(aes(x = time, color = type)) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Objects [satellites]",
      color = "Object type"
    ) +
    scale_color_manual(
      name = "Object type",
      breaks = c("Derelict", "Constellation", "Fringe"),
      labels = c("Derelict", "Constellation", "Fringe"),
      values = c("Derelict" = color_set[1], "Constellation" = color_set[2], "Fringe" = color_set[3])
    )
  }

  ymax = df %>% filter(variable==variable_name) %>% ungroup() %>% select(values) %>% max()
  ymin = df %>% filter(variable==variable_name) %>% ungroup() %>% select(values) %>% min()

  plot_out = plot_base +
    geom_line( aes(y = values), linewidth = 1) +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 15)
    ) +
    theme_bw() +
    ylim(c(ymin, ymax))

  if(leftmost==FALSE) { plot_out = plot_out + labs(y="") }

  if(variable_name=="expected_max_loss") {
    plot_out = plot_out +
    scale_y_continuous( # Set label and then set ymin and ymax to ymin and ymax
      labels = label_number(scale = 1e-6, big.mark = ","),
      limits = c(ymin, ymax)
      )
  }

  plot_out
}

# Non-normalized aggregate SSR index
ssr_left = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_1, 20, scenario_label_1, "ssr_index", TRUE)
ssr_right = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_2, 20, scenario_label_2, "ssr_index", FALSE)

combined_ssr_plot = ssr_left + ssr_right + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
combined_ssr_plot = combined_ssr_plot & theme(legend.position = "bottom")

ggsave(
  paste0("scenarios/", out_folder, "/ssr-comparison.png"), 
  bg = "white",
  combined_ssr_plot, width = 8, height = 6, dpi = 320)

# Normalized aggregate SSR index
ssr_left_norm = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_1, 20, scenario_label_1, "normalized_ssr_index", TRUE)
ssr_right_norm = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_2, 20, scenario_label_2, "normalized_ssr_index", FALSE)

combined_ssr_plot_norm = ssr_left_norm + ssr_right_norm + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
combined_ssr_plot_norm = combined_ssr_plot_norm & theme(legend.position = "bottom")

ggsave(
  paste0("scenarios/", out_folder, "/normalized-ssr-comparison.png"), 
  bg = "white",
  combined_ssr_plot_norm, width = 8, height = 6, dpi = 320)

# Expected maximum losses
emaxloss_left = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_1, 20, scenario_label_1, "expected_max_loss", TRUE)
emaxloss_right = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_2, 20, scenario_label_2, "expected_max_loss", FALSE)

combined_emaxloss_plot = emaxloss_left + emaxloss_right + plot_annotation(tag_levels = 'A', title="")

ggsave(
  paste0("scenarios/", out_folder, "/expected-maximum-welfare-comparison.png"), 
  bg = "white",
  combined_emaxloss_plot, width = 8, height = 6, dpi = 320)

# Total number of sats
total_sats_left = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_1, 20, scenario_label_1, "total_objects", TRUE)
total_sats_right = comp_plot_fn(aggregate_metrics_comparison_long, scenario_label_2, 20, scenario_label_2, "total_objects", FALSE)

combined_total_sats_plot = total_sats_left + total_sats_right + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
combined_total_sats_plot = combined_total_sats_plot & theme(legend.position = "bottom")

ggsave(
  paste0("scenarios/", out_folder, "/total_sats-comparison.png"), 
  bg = "white",
  combined_total_sats_plot, width = 8, height = 6, dpi = 320)


# # Plot total collision probability for each type over time. Replace "type" with "Type" in the legend, and replace labels as D=Derelict, S=Slotted, Su=Unslotted.
# ## Generate plot
# plot_total_collision_probability_D = aggregate_metrics_comparison_long %>%
#   filter(type == "D") %>%
#   ggplot(aes(x = time, color = label)) +
#   geom_line( aes(y = total_collision_probability), linewidth = 1) +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.text.x = element_text(size = 15),
#     axis.text.y = element_text(size = 15),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 12),
#     plot.title = element_text(size = 15)
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Difference in total collision probability over time [Derelicts]",
#     x = "Time",
#     y = "[pp]"
#   )

# ggsave(
#   paste0("scenarios/", stem_1, "/", input_name_1,"__total-collision-probability-D.png"), 
#   bg = "white",
#   plot_total_collision_probability_D, width = 8, height = 6, dpi = 320)

# plot_total_collision_probability_Su = aggregate_metrics_comparison_long %>%
#   filter(type == "Su") %>%
#   ggplot(aes(x = time, color = label)) +
#   geom_line( aes(y = total_collision_probability), linewidth = 1) +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.text.x = element_text(size = 15),
#     axis.text.y = element_text(size = 15),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 12),
#     plot.title = element_text(size = 15)
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Difference in total collision probability over time [Fringe]",
#     x = "Time",
#     y = "[pp]"
#   )

# ggsave(
#   paste0("scenarios/", stem_1, "/", input_name_1,"__total-collision-probability-Su.png"), 
#   bg = "white",
#   plot_total_collision_probability_Su, width = 8, height = 6, dpi = 320)


# plot_total_collision_probability_S = aggregate_metrics_comparison_long %>%
#   filter(type == "S") %>%
#   ggplot(aes(x = time, color = label)) +
#   geom_line( aes(y = total_collision_probability), linewidth = 1) +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.text.x = element_text(size = 15),
#     axis.text.y = element_text(size = 15),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 12),
#     plot.title = element_text(size = 15)
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Difference in total collision probability over time [Constellation]",
#     x = "Time",
#     y = "[pp]"
#   )

# ggsave(
#   paste0("scenarios/", stem_1, "/", input_name_1,"__total-collision-probability-S.png"), 
#   bg = "white",
#   plot_total_collision_probability_S, width = 8, height = 6, dpi = 320)

