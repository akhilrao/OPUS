#!/bin/bash

### 1. INITIAL DEFINITIONS

## Don't need to change, defined here for reducing typo errors later
launch_pattern_type_equilibrium="equilibrium"
launch_pattern_type_feedback="sat_feedback"
model_type="MOCAT"

# Define the length of the horizon over which to solve (years)
model_horizon=35

# Define the R script(s) to run
R_SCRIPT="analytics.R" 

### 2. SCENARIO DEFINITIONS
# Select the scenarios to process. Store them in an array. Directly use human-readable names here, this script only runs R files.
scenario_names=(
                "staying-coral-stew-Soybean"
                )

### 3.  EXECUTION
# Loop over scenarios with the scenario_file inputs and run the analytics script.

for scenario_name in "${scenario_names[@]}"
do
    # Launch behavior 1: Equilibrium, MOCAT or GMPHD
    Rscript "analytics.R" $scenario_name $model_type $launch_pattern_type_equilibrium $model_horizon

    # Launch behavior 2: Feedback, MOCAT only
    Rscript "analytics.R" $scenario_name $model_type $launch_pattern_type_feedback $model_horizon

    echo "Scenario ~~ $scenario_name ~~ processing complete!"
done