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
library(scales)
source("analytics-functions.R")

# Read in command line arguments: stem, model_type, and launch_pattern_type, and model_horizon
args = commandArgs(trailingOnly = TRUE)
stem = args[1]
model_type = args[2]
launch_pattern_type = args[3]
model_horizon = args[4]

# Script version
# stem = "destroying-termite-verbalize-Zinnia"
# model_type = "GMPHD"
# launch_pattern_type = "equilibrium"
# model_horizon = "10"

# Construct input filename
input_name = paste0(stem, "-", model_type, "-", launch_pattern_type, "-", model_horizon, "yrs") #Input filename

# Read the files and compute altitude_bin_midpoints
altitude_bin_midpoints = read_and_compute_midpoints("x0_TLE/Counts_DEBRIS_bins_35.csv")

# Read the output CSV
output_data = read_csv(paste0("scenarios/", stem, "/", input_name, "__out-data.csv"))

# Read the collision parameters CSV
collision_parameters = cbind(
  altitude_bin_midpoints,
  read_csv(paste0("scenarios/", stem, "/", stem, "--collision-parameters.csv"))
  )

# Read the fragmentation parameters CSV
fragmentation_parameters = read_csv(paste0("scenarios/", stem, "/fragmentation-parameters.csv"))

# Read the parameters CSV
parameters = read_csv(paste0("scenarios/", stem, "/", stem, "--parameters.csv"))
## Calculate costs and revenues for individual satellite
### costs
cost = parameters %>% filter(parameter_type=="econ", parameter_name=="cost") %>% select(c(3:42)) %>% t() %>% as.numeric()
cost
stationkeeping_cost = parameters %>% filter(parameter_type=="econ", parameter_name=="stationkeeping_cost") %>% select(-(c(1,2))) %>% t() %>% as.data.frame() %>% drop_na() %>% unlist() %>% as.numeric()
stationkeeping_cost
lifetime_loss_cost = parameters %>% filter(parameter_type=="econ", parameter_name=="lifetime_loss_cost") %>% select(-(c(1,2))) %>% t() %>% as.data.frame() %>% drop_na() %>% unlist() %>% as.numeric()
lifetime_loss_cost
### revenues
revenues_intercept = parameters %>% filter(parameter_type=="econ", parameter_name=="intercept") %>% select(3) %>% t() %>% as.numeric()
revenues_slope = parameters %>% filter(parameter_type=="econ", parameter_name=="coef") %>% select(3) %>% t() %>% as.numeric()
#### To get the fringe satellite patterns, select the "Su_" columns from the output data, then compute sums over each column
fringe_satellite_patterns = output_data %>% select(starts_with("Su_")) %>% colSums()
revenues = revenues_intercept - revenues_slope * fringe_satellite_patterns
### Accounting rate of return ignoring collision risk. Give column names that are "gross_ror_" followed by the index number of the column.
accounting_profits = matrix(revenues, nrow=length(cost), ncol=length(revenues), byrow=TRUE)/matrix(cost, nrow=length(cost), ncol=length(revenues), byrow=FALSE)
colnames(accounting_profits) = paste0("gross.ror_", 1:length(revenues))
## Extract discount rate
discount_rate = parameters %>% filter(parameter_type=="econ", parameter_name=="discountRate") %>% select(3) %>% unlist() %>% as.numeric()
## Calculate NPV of satellite assuming no collisions. This is an upper bound on the private WTP to replace a satellite that was just launched
pv_revenues_time = discounted_revenues(revenues, discount_rate)
collision_risk_free_npv = matrix(pv_revenues_time, nrow=length(cost), ncol=length(pv_revenues_time), byrow=TRUE) - matrix(cost, nrow=length(cost), ncol=length(pv_revenues_time), byrow=FALSE)
colnames(collision_risk_free_npv) = paste0("risk.free.npv_", 1:length(pv_revenues_time))

# Add altitude_bin_midpoints to output_data
output_data_aug = data.frame(altitude_bin_midpoints, cost_0=cost, skCost_0=stationkeeping_cost, llCost_0=lifetime_loss_cost, accounting_profits, collision_risk_free_npv, output_data)

# Reshape data for ggplot
output_data_melted = output_data_aug %>%
  pivot_longer(cols = -altitude_bin_midpoints, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("type", "time"), sep = "_") %>%
  mutate(time = as.numeric(time))

if(model_type=="GMPHD") {
  plot_S = heatmap_plot(output_data_melted, "S", "Constellation sats")
  plot_Su = heatmap_plot(output_data_melted, "Su", "Fringe satellites")
  plot_Dlt = heatmap_plot(output_data_melted, "Dlt", "Large trackable debris")
  plot_Dst = heatmap_plot(output_data_melted, "Dst", "Small trackable debris")
  plot_Dlnt = heatmap_plot(output_data_melted, "Dlnt", "Lethal nontrackables")

  final_plot = (plot_Dlt | plot_Dst | plot_Dlnt) + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
  final_plot = final_plot & theme(legend.position = "bottom")

  ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots.png"), final_plot, width = 16, height = 12, dpi = 320, bg="white")

  # save separate heatmap plots
  plot_S_solo = heatmap_plot(output_data_melted, "S", "Constellation sats", bar=1)
  plot_Su_solo = heatmap_plot(output_data_melted, "Su", "Fringe satellites", bar=1)
  plot_Dlt_solo = heatmap_plot(output_data_melted, "Dlt", "Large trackable debris", bar=1)
  plot_Dst_solo = heatmap_plot(output_data_melted, "Dst", "Small trackable debris", bar=1)
  plot_Dlnt_solo = heatmap_plot(output_data_melted, "Dlnt", "Lethal nontrackables", bar=1)

  ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-S.png"),  plot_S_solo, width = 5, height = 3, dpi = 320, bg="white")
  ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-Su.png"),  plot_Su_solo, width = 5, height = 3, dpi = 320, bg="white")
  ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-Dlt.png"),  plot_Dlt_solo, width = 5, height = 3, dpi = 320, bg="white")
  ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-Dst.png"),  plot_Dst_solo, width = 5, height = 3, dpi = 320, bg="white")
  ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-Dlnt.png"),  plot_Dlnt_solo, width = 5, height = 3, dpi = 320, bg="white")

}

# Loop over time periods within types to fill in new collision_probability and ecob variables, initialized as NA.
output_data_melted$collision_probability = NA
output_data_melted$ecob = NA
for(i in 1:length(output_data_melted$time)) {
  df_t = output_data_melted[output_data_melted$time==output_data_melted$time[i],]

  # Constellation
  output_data_melted$collision_probability[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="S"] = collisionProbability(df_t, collision_parameters, sat_type="S", rate_type="negative_exponential")
  output_data_melted$ecob[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="S"] = ecob_index(df_t, collision_parameters, fragmentation_parameters, type="S")

  # Fringe
  output_data_melted$collision_probability[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="Su"] = collisionProbability(df_t, collision_parameters, sat_type="Su", rate_type="negative_exponential")
  output_data_melted$ecob[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="Su"] = ecob_index(df_t, collision_parameters, fragmentation_parameters, type="Su")

  # Derelict
  output_data_melted$collision_probability[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="D"] = collisionProbability(df_t, collision_parameters, sat_type="D", rate_type="negative_exponential")
  output_data_melted$ecob[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="D"] = ecob_index(df_t, collision_parameters, fragmentation_parameters, type="D")
}

# Attach expected loss as new column
output_data_melted_2 = output_data_melted %>% filter(type=="S" | type=="Su" | type=="D")
output_data_melted_2$risk.free.npv = NA
output_data_melted_2$risk.free.npv[output_data_melted_2$type=="Su"] = output_data_melted %>% filter(type=="risk.free.npv") %>% select(value) %>% unlist()

# Compute total collision_probability for each type and each across all altitude_bin_midpoints
aggregate_metrics_paths = output_data_melted_2 %>%
  filter(type=="S" | type=="Su" | type=="D") %>% 
  group_by(time) %>%
  mutate(
    ecob = sum(ecob, na.rm = TRUE)
  ) %>%
  group_by(type, time) %>%
  summarise(
    total_objects = sum(value, na.rm=TRUE),
    total_collision_probability = sum(collision_probability, na.rm = TRUE),
    mean_pc = mean(collision_probability, na.rm = TRUE),
    ecob = sum(ecob, na.rm = TRUE),
    expected_max_loss = sum(value*pmax((1-collision_probability)*risk.free.npv,0),na.rm=TRUE) # take pmax here to reflect that firms won't launch satellites that won't be profitable even without collisions
    ) %>%
  mutate(
    ecob = ecob/ecob[1],
    ssr_index = total_collision_probability*ecob,
    normalized_ssr_index = ssr_index/ssr_index[1]
  )

aggregate_metrics_paths %>% print(n=Inf)

write_csv(aggregate_metrics_paths, paste0("scenarios/", stem, "/", input_name, "__aggregate-metrics-paths.csv"))

#####
# Generate plots
#####

# Plot total collision probability for each type over time. Replace "type" with "Type" in the legend, and replace labels as D=Derelict, S=Constellation, Su=Fringe.
## Extract the first three colors from the Set2 palette
color_set <- brewer.pal(3, "Set1")
## Generate plot
plot_total_collision_probability = aggregate_metrics_paths %>%
  ggplot(aes(x = time, y = total_collision_probability, color = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    name = "Type",
    breaks = c("D", "S", "Su"),
    labels = c("Derelict", "Constellation", "Fringe"),
    values = c("D" = color_set[1], "S" = color_set[2], "Su" = color_set[3])
  ) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "",
    x = "Time",
    y = "Annual expected collisions"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__total-collision-probability.png"), 
  bg = "white",
  plot_total_collision_probability, width = 8, height = 6, dpi = 320)

plot_ecob = aggregate_metrics_paths %>% filter(type=="S") %>%
  ggplot(aes(x = time, y = ecob)) +
  geom_line(linewidth = 1) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "",
    x = "Time",
    y = "ECOB Index"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__ecob.png"), 
  bg = "white",
  plot_ecob, width = 8, height = 6, dpi = 320)


plot_ssr_index = aggregate_metrics_paths %>%
  ggplot(aes(x = time, y = ssr_index, color = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    name = "Type",
    breaks = c("D", "S", "Su"),
    # labels = c("Derelict", "Constellation", "Fringe"),
    values = c("D" = color_set[1], "S" = color_set[2], "Su" = color_set[3])
  ) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "",
    x = "Time",
    y = "SSR Index"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__ssr-index.png"), 
  bg = "white",
  plot_ssr_index, width = 8, height = 6, dpi = 320)


# Make and combine heatmap plots
plot_S = heatmap_plot(output_data_melted, "S", "Constellation sats")
plot_Su = heatmap_plot(output_data_melted, "Su", "Fringe satellites")
plot_N = heatmap_plot(output_data_melted, "N", "Debris")
plot_D = heatmap_plot(output_data_melted, "D", "Derelict satellites")
plot_launch_slotted = heatmap_plot(output_data_melted, "lConstellation", "Constellation launches", leftmost=TRUE)
plot_launch_unslotted = heatmap_plot(output_data_melted, "lFringe", "Fringe launches", leftmost=TRUE)

final_plot = (plot_launch_slotted + plot_S + plot_N) / (plot_launch_unslotted + plot_Su + plot_D) + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
final_plot = final_plot & theme(legend.position = "bottom")

ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots.png"), final_plot, width = 16, height = 12, dpi = 320, bg="white")

final_plot_Su_D = (plot_Su + plot_D) + plot_annotation(tag_levels = 'A', title="") + plot_layout(guides = 'collect')
final_plot_Su_D = final_plot_Su_D & theme(legend.position = "bottom")

ggsave(paste0("scenarios/", stem, "/Su-D-", input_name,"__heatmap-time-plots.png"), final_plot_Su_D, width = 10*0.75, height = 6*0.75, dpi = 320, bg="white")

# save separate heatmap plots for derelicts, fringe satellites, constellation satellites, fragments over time
plot_Su_solo = heatmap_plot(output_data_melted, "Su", "Fringe sats", barwidth=1)
plot_D_solo = heatmap_plot(output_data_melted, "D", "Derelict sats", barwidth=1)
plot_S_solo = heatmap_plot(output_data_melted, "S", "Constellation sats", barwidth=1)
plot_N_solo = heatmap_plot(output_data_melted, "N", "Debris", barwidth=1)

ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-D.png"), plot_D_solo, width = 5, height = 3, dpi = 320, bg="white")
ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-Su.png"),  plot_Su_solo, width = 5, height = 3, dpi = 320, bg="white")
ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-S.png"),  plot_S_solo, width = 5, height = 3, dpi = 320, bg="white")
ggsave(paste0("scenarios/", stem, "/", input_name,"__heatmap-time-plots-N.png"),  plot_N_solo, width = 5, height = 3, dpi = 320, bg="white")

# Plot scatterplot of cost function and components (stationkeeping, lifetime loss) over locations
plot_cost_scatter  = output_data_melted %>%
  filter(type=="cost") %>%
  ggplot(aes(x = altitude_bin_midpoints, y = value)) +
  geom_point() +
  scale_y_continuous(labels = label_number(scale = 1e-6, big.mark = ",")) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "",
    x = "Altitude [km]",
    y = "Cost [M$/sat]"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__scatter-cost.png"), plot_cost_scatter, 
  bg="white",
  width = 8, 
  height = 6, 
  dpi = 320)

plot_cost_line = output_data_melted %>%
  filter(type=="cost") %>%
  ggplot(aes(x = altitude_bin_midpoints, y = value)) +
  geom_line(linewidth=1) +
  scale_y_continuous(labels = label_number(scale = 1e-6, big.mark = ",")) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  theme_bw() +
  labs(
    title = "",
    x = "Altitude [km]",
    y = "Cost [M$/sat]"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__line-cost.png"), plot_cost_line, 
  bg="white",
  width = 6, 
  height = 6, 
  dpi = 320)

plot_stationkeeping_cost_scatter = output_data_melted %>%
  filter(type=="skCost") %>%
  ggplot(aes(x = altitude_bin_midpoints, y = value)) +
  geom_point() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "Stationkeeping cost over altitude",
    x = "Altitude",
    y = "Stationkeeping cost"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__scatter-stationkeeping-cost.png"), plot_stationkeeping_cost_scatter, 
  bg="white",
  width = 8, 
  height = 6, 
  dpi = 320)

plot_lifetime_loss_cost_scatter = output_data_melted %>%
  filter(type=="llCost") %>%
  ggplot(aes(x = altitude_bin_midpoints, y = value)) +
  geom_point() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "Lifetime loss cost over altitude",
    x = "Altitude",
    y = "Lifetime loss cost"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__scatter-lifetime-loss-cost.png"), plot_lifetime_loss_cost_scatter, 
  bg="white",
  width = 8, 
  height = 6, 
  dpi = 320)

# Plot scatterplot of accounting gross rate of return over locations in the first column, first in period 0 then in all periods
plot_accounting_profits_scatter_0 = output_data_melted %>%
  filter(type=="gross.ror", time==1) %>%
  ggplot(aes(x = altitude_bin_midpoints, y = value*100)) +
  geom_point() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  labs(
    title = "Accounting rate of return over altitude",
    x = "Altitude",
    y = "[%]"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__scatter-accounting-profits-0.png"), plot_accounting_profits_scatter_0, 
  bg="white",
  width = 8, 
  height = 6, 
  dpi = 320)

plot_accounting_profits_scatter_all = output_data_melted %>%
  filter(type=="gross.ror") %>%
  ggplot(aes(x = altitude_bin_midpoints, y = value*100, color = time)) +
  geom_point() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  ) +
  theme_minimal() +
  scale_color_distiller(palette = "YlOrRd") +
  labs(
    title = "Accounting rate of return over altitude",
    x = "Altitude",
    y = "[%]"
  )

ggsave(
  paste0("scenarios/", stem, "/", input_name,"__scatter-accounting-profits-all.png"), plot_accounting_profits_scatter_all, 
  bg="white",
  width = 8, 
  height = 6, 
  dpi = 320)

# Write out summary statistics of interest in a csv file, columns are labels and rows are entries:
# Location of maximum value for each object type and maximum value
# Average launch rate of slotted and unslotted satellites across all bins and all time
# Location of maximum cumulative annual growth rate of debris and derelict satellites, maximum cumulative annual growth rate across bins, and harmonic mean cumulative annual growth rate across all bins and all time

summary_stats = output_data_melted %>%
  group_by(type, altitude_bin_midpoints) %>%
  mutate(cagr = ((value[n()] / value[1]) ^ (1 / (time[n()] - time[1])) - 1)*100 ) %>%
  ungroup() %>% group_by(type) %>%
  summarise(
    max_value = round(max(value, na.rm = TRUE),2),
    time_of_max = time[which.max(value)],
    location_of_max = altitude_bin_midpoints[which.max(value)],
    avg_launch_rate = round(mean(value[grepl("l", type)], na.rm = TRUE),2),
    max_cagr = round(max(cagr, na.rm = TRUE),2),
    time_of_max_cagr = time[which.max(cagr)],
    location_of_max_cagr = altitude_bin_midpoints[which.max(cagr)],
    harmonic_mean_cagr = round(1 / mean(1 / cagr, na.rm = TRUE),2)
    )

write_csv(summary_stats, paste0("scenarios/", stem, "/", input_name, "__summary-stats.csv"))