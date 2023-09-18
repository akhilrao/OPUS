#####
# R functions to generate figures and calculate statistics/figures of merit about the simulation. 
# "Slotted" and "Unslotted" are used as synonyms for "Constellation" and "Fringe" respectively.
#####

### HELPER FUNCTION BLOCK ###

safe_divide <- function(a, b) {
  # Vectorized function to divide two numbers, but return 0 if a + b is 0
  # Input: two vectors of numbers
  # Output: vector of quotients of the two numbers, or 0 where a + b is 0
  result <- ifelse((a + b) == 0, 0, a / (a + b))
  return(result)
}

read_and_compute_midpoints <- function(filename) {
  # Function to read CSV and compute midpoints
  # Input: filename of CSV
  # Output: vector of altitude bin midpoints

  data = read_csv(filename)
  
  # Extract the numerical values from the range column
  ranges = str_extract_all(data$range, "\\d+\\.\\d+") %>%
    unlist() %>%
    as.numeric() %>%
    matrix(ncol = 2, byrow = TRUE)
  
  # Compute the midpoints of each range, and subtract the constant
  altitude_bin_midpoints = rowMeans(ranges) - 6378.1366
  
  return(altitude_bin_midpoints)
}

# Collision probability function definitions for objects of type S, Su, N. This is contemporaneous collision probability, not next-period. 
# Each collision probability is constructed according to two formulas: 
# negative exponential: 1 - exp(-collision_rate)
# linear fraction: collision_rate/number_of_satellites, e.g. for slotted "slotted_collision_rate/number_of_slotted_satellites"
# Function to compute collision probability for any type of satellite
collisionProbability <- function(df, collision_parameters, sat_type="S", rate_type="negative_exponential") {
# Input: collision_parameters data frame, output_data_melted data frame
# Output: vector of collision probabilities for slotted satellites

  # Define variables for each object type
  S_value  = df$value[df$type == "S"]
  N_value  = df$value[df$type == "N"]
  D_value  = df$value[df$type == "D"]
  Su_value = df$value[df$type == "Su"]

  # Define collision parameters for each interaction type
  SN = abs(collision_parameters$SlottedDebris)
  SD = abs(collision_parameters$SlottedDerelict)
  SSu = abs(collision_parameters$SlottedUnslotted)
  SS = abs(collision_parameters$SlottedSlotted)
  SuN = abs(collision_parameters$UnslottedDebris)
  SuD = abs(collision_parameters$UnslottedDerelict)
  SuS = abs(collision_parameters$UnslottedSlotted)
  SuSu = abs(collision_parameters$UnslottedUnslotted)
  DN = abs(collision_parameters$DerelictDebris)
  DD = abs(collision_parameters$DerelictDerelict)
  DS = abs(collision_parameters$DerelictSlotted)
  DSu = abs(collision_parameters$DerelictUnslotted)

  # Compute density of interactions with objects large enough to trigger a collision. These are the -\rho \Delta V A \Delta t terms in equation (4) of Letizia et al 2017.
  if(sat_type=="S"){
      collision_rate = SN  * S_value * N_value +
        SD * S_value * D_value +
        SSu * S_value * Su_value * safe_divide(S_value, Su_value) +
        SS * S_value ^ 2                                
    }
  if(sat_type=="Su"){
      collision_rate = SuN  * Su_value * N_value +
        SuD * Su_value * D_value +
        SuS * Su_value * S_value * safe_divide(Su_value, S_value) +
        SuSu * Su_value ^ 2
    }
  if(sat_type=="D"){
      collision_rate = DN  * D_value * N_value +
        DS * D_value * S_value +
        DSu * D_value * Su_value +
        DD * D_value ^ 2
    }

  # Compute negative exponential collision probability if switched on
  if(rate_type=="negative_exponential"){
    collision_probability = 1 - exp(-collision_rate)
  }

  return(collision_probability)
}

# Function to compute ECOB index for Su objects. ECOB formula is: \sum_{j=1}^{N} w_j prob_j, where j=1,,,N represent the different types of objects being considered (S, Su, D). The algorithm is as follows:
# 1. Select a target set composed by N objects
# 2. For each object under evaluation,
#     a. Compute its collision probability (use pre-computed collision probabilities from the collisionProbability() function)
#     b. Simulate its fragmentation (use pre-computed fragmentation parameters from the K0.csv object)
#     c. Evaluate the change in Pc for each object in N due to the simulated fragmentation
#     d. Sum the contribution for each object in N (using cumulative Pc facing that object post-fragmentation and weighting factor)
# 3. Obtain object-specific ECOB numbers that measure the risk (prob * consequence) of that object breaking up
ecob_index <- function(df, collision_parameters, fragmentation_parameters, type="S") {
# Input: output_data_melted data frame, collision_parameters data frame, fragmentation_parameters data frame

  # Compute cross-sectional area ratios. S, Su, D type objects have areas of 1.741, N type objects have areas of 0.02. For each altitude compute the total amount of cross-sectional area as number_of_objects*area (summed over all object types). Then group by altitude and compute the share of cross-sectional area at each altitude.
  weights = df %>% 
    filter(type=="Su" | type=="S" | type=="D") %>%
    group_by(altitude_bin_midpoints) %>%
    summarise(total_area = sum(value * ifelse(type=="Su" | type=="S" | type=="D", 1.741, 0.02))) %>%
    mutate(weight = total_area/sum(total_area)) %>%
    select(altitude_bin_midpoints, weight)

  # Compute Pc
  pc = collisionProbability(df, collision_parameters, sat_type=type, rate_type="negative_exponential")

  # Simulate fragmentation for Su type objects
  new_fragments = max(fragmentation_parameters[,type]) # Since the question here is about catastrophic fragmentations, take the largest value in the Su column -- under MOCAT4S this is common for all fragmentation pairs except with small objects
  df_updated = df
  df_updated$value[df$type == "N"] = df$value[df$type == "N"] + new_fragments*(df$value[df$type == type]>1) # Compute the updated number of small fragments objects after fragmentation
  df_updated$value[df$type == type] = pmax(df$value[df$type == type] - 1,0) # Compute the updated number of Su objects after fragmentation

  # Update Pc
  pc_updated = collisionProbability(df_updated, collision_parameters, sat_type=type, rate_type="negative_exponential")

  # Compute ECOB
  ecob = weights$weight * pc_updated

  return(ecob)
}

# Function to compute economic welfare. 

# Function to prepare data for comparison
prepare_data_for_comparison <- function(stem, input_name, altitude_bin_midpoints) {
# Input: stem, input_name, altitude_bin_midpoints
# Output: list of data frames

  # Read the output CSV for scenario
  output_data = read_csv(paste0("scenarios/", stem, "/", input_name, "__out-data.csv"))

  # Read the collision parameters CSV
  collision_parameters = cbind(
    altitude_bin_midpoints,
    read_csv(paste0("scenarios/", stem, "/", stem, "--collision-parameters.csv"))
    )

  # Read the parameters CSV
  parameters = read_csv(paste0("scenarios/", stem, "/", stem, "--parameters.csv"))
  ## Calculate costs and revenues for individual satellite
  ### costs
  cost = parameters %>% filter(parameter_type=="econ", parameter_name=="cost") %>% select(c(3:42)) %>% t() %>% as.numeric()
  cost
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

  # Add altitude_bin_midpoints to output_data
  output_data_aug = data.frame(altitude_bin_midpoints, accounting_profits, output_data)

  # Reshape data for ggplot
  output_data_melted = output_data_aug %>%
    pivot_longer(cols = -altitude_bin_midpoints, names_to = "variable", values_to = "value") %>%
    separate(variable, into = c("type", "time"), sep = "_") %>%
    mutate(time = as.numeric(time))

  # Loop over time periods within types to fill in a new variable labeled collision_probability. This variable is initialized as NA. Then, for each time period, compute the collision probability for slotted satellites (type=="S") and assign it to the appropriate rows.
  output_data_melted$collision_probability = NA
  for(i in 1:length(output_data_melted$time)) {
    # Slotted
    output_data_melted$collision_probability[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="S"] = collisionProbability(output_data_melted[output_data_melted$time==output_data_melted$time[i],], collision_parameters, sat_type="S", rate_type="negative_exponential")
    # Unslotted
    output_data_melted$collision_probability[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="Su"] = collisionProbability(output_data_melted[output_data_melted$time==output_data_melted$time[i],], collision_parameters, sat_type="Su", rate_type="negative_exponential")
    # Derelict
    output_data_melted$collision_probability[output_data_melted$time==output_data_melted$time[i] & output_data_melted$type=="D"] = collisionProbability(output_data_melted[output_data_melted$time==output_data_melted$time[i],], collision_parameters, sat_type="D", rate_type="negative_exponential")
  }

  # Compute total collision_probability for each type and each across all altitude_bin_midpoints
  aggregate_metrics_paths = output_data_melted %>%
    filter(type=="S" | type=="Su" | type=="D") %>%
    group_by(type, time) %>%
    summarise(total_collision_probability = sum(collision_probability, na.rm = TRUE))

  return(list(output_data=output_data, collision_parameters=collision_parameters, parameters=parameters, cost=cost, fringe_satellite_patterns=fringe_satellite_patterns, revenues=revenues, accounting_profits=accounting_profits, discount_rate=discount_rate, output_data_aug=output_data_aug, output_data_melted=output_data_melted, aggregate_metrics_paths=aggregate_metrics_paths))
}

heatmap_plot <- function(data, obj_type, type_name, leftmost=FALSE, barwidth=10, model_type="MOCAT") {
  # Function to generate heatmap for each object type. Should only have a y axis title if "leftmost==TRUE"
  # Input: data frame, object type, type name, whether or not it is the leftmost plot, and the width of the colorbar
  # Output: ggplot object

    min_value = min(data$value, na.rm = TRUE)
    max_value = max(data$value, na.rm = TRUE)

  if(model_type=="MOCAT"){
    # Check if the "type_name" value contains "Fringe"; if it does, assign "palette_val" to "Blues". If it does not, check if it contains "Constellation"; if it does, assign it to "Reds". Otherwise, assign it to "Greens". In all cases, assign "min_value" and "max_value" to the minimum and maximum values of "value" for the rows where "type" contains the key used in "type_name"; e.g. if type_name=="Fringe" then assign "min_value" based on the rows where type contains "Su" or "lUnslotted". Finally, assign "legend_name" to the appropriate value based on "type_name"
    if(grepl("Fringe", type_name)) {
      palette_val = "Blues"
      min_value = min(data$value[data$type == "Su" | data$type == "lUnslotted"], na.rm = TRUE)
      max_value = max(data$value[data$type == "Su" | data$type == "lUnslotted"], na.rm = TRUE)
      legend_name = "Fringe"
    } else if(grepl("Constellation", type_name)) {
      palette_val = "YlGn"
      min_value = min(data$value[data$type == "S" | data$type == "lSlotted"], na.rm = TRUE)
      max_value = max(data$value[data$type == "S" | data$type == "lSlotted"], na.rm = TRUE)
      legend_name = "Constellation"
    } else {
      palette_val = "OrRd"
      min_value = min(data$value[data$type == "N" | data$type == "D"], na.rm = TRUE)
      max_value = max(data$value[data$type == "N" | data$type == "D"], na.rm = TRUE)
      legend_name = "Junk"
    }
  }

  if(model_type=="GMPHD"){
    # Check if the "type_name" value contains "Fringe"; if it does, assign "palette_val" to "Blues". If it does not, check if it contains "Constellation"; if it does, assign it to "Reds". Otherwise, assign it to "Greens". In all cases, assign "min_value" and "max_value" to the minimum and maximum values of "value" for the rows where "type" contains the key used in "type_name"; e.g. if type_name=="Fringe" then assign "min_value" based on the rows where type contains "Su" or "lUnslotted". Finally, assign "legend_name" to the appropriate value based on "type_name"
    if(grepl("Fringe", type_name)) {
      palette_val = "Blues"
      min_value = min(data$value[data$type == "Su" | data$type == "lUnslotted"], na.rm = TRUE)
      max_value = max(data$value[data$type == "Su" | data$type == "lUnslotted"], na.rm = TRUE)
      legend_name = "Fringe"
    } else if(grepl("Constellation", type_name)) {
      palette_val = "YlGn"
      min_value = min(data$value[data$type == "S" | data$type == "lSlotted"], na.rm = TRUE)
      max_value = max(data$value[data$type == "S" | data$type == "lSlotted"], na.rm = TRUE)
      legend_name = "Constellation"
    } else {
      palette_val = "OrRd"
      min_value = min(data$value[data$type == "Dlt" | data$type == "Dst" | data$type == "Dlnt"], na.rm = TRUE)
      max_value = max(data$value[data$type == "Dlt" | data$type == "Dst" | data$type == "Dlnt"], na.rm = TRUE)
      legend_name = "Junk"
    }
  }

  plot_out = data %>%
    filter(type == obj_type) %>%
    ggplot(aes(x = time, y = altitude_bin_midpoints, fill = value)) +
    geom_tile() +
    scale_fill_distiller(
      palette = palette_val, 
      limits = c(min_value, max_value), 
      guide = guide_colorbar(barwidth = barwidth), 
      direction=1,
      name=legend_name
      ) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 15)) +
    theme_bw()

  if(leftmost==FALSE) {
  plot_out = plot_out +
    labs(title = type_name, x="Time", y="")
  }
  if(leftmost==TRUE) {
  plot_out = plot_out +
    labs(title = type_name, x="Time", y="Altitude")
  }

  plot_out
}


discounted_revenues <- function(revenues, discount_rate=0.05) { 
  # Function to calculate present value of revenues for a satellite over 5-year lifespan
  
  pv_revenues_mat = data.frame(revenues_0=revenues) %>%
    mutate(
      revenues_1 = data.table::shift(revenues, n=1, type="lead")*(1/(1+discount_rate)),
      revenues_2 = data.table::shift(revenues, n=2, type="lead")*(1/(1+discount_rate))^2,
      revenues_3 = data.table::shift(revenues, n=3, type="lead")*(1/(1+discount_rate))^3,
      revenues_4 = data.table::shift(revenues, n=4, type="lead")*(1/(1+discount_rate))^4
    )
    for(r in 1:nrow(pv_revenues_mat)) {
      pv_revenues_mat[r,"revenues_1"] = ifelse(is.na(pv_revenues_mat$revenues_1[r]), pv_revenues_mat$revenues_0[r]*(1/(1+discount_rate)), pv_revenues_mat$revenues_1[r])
      pv_revenues_mat[r,"revenues_2"] = ifelse(is.na(pv_revenues_mat$revenues_2[r]), pv_revenues_mat$revenues_1[r]*(1/(1+discount_rate)), pv_revenues_mat$revenues_2[r])
      pv_revenues_mat[r,"revenues_3"] = ifelse(is.na(pv_revenues_mat$revenues_3[r]), pv_revenues_mat$revenues_2[r]*(1/(1+discount_rate)), pv_revenues_mat$revenues_3[r])
      pv_revenues_mat[r,"revenues_4"] = ifelse(is.na(pv_revenues_mat$revenues_4[r]), pv_revenues_mat$revenues_3[r]*(1/(1+discount_rate)), pv_revenues_mat$revenues_4[r])
    }
    pv_revenues = rowSums(pv_revenues_mat)
    return(pv_revenues)
}