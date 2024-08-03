#!/bin/bash

### 1. INITIAL DEFINITIONS

## Don't need to change, defined here for reducing typo errors later
launch_pattern_type_equilibrium="equilibrium"
launch_pattern_type_feedback="sat_feedback"
model_type="MOCAT"

# Define the R script(s) to run
R_SCRIPT="compare-two-scenarios.R" 

## Define a folder that the comparisons will be saved under
mkdir scenarios/comparison-benchmark-eqm-satfeedback

### 2. SCENARIO DEFINITIONS
# Define the length of the horizon over which to solve (years)
model_horizon=35

## Specify first scenario for comparison
scenario_name_1="generous-turtle-tastes-Hibiscus"
                
## Specify second scenario for comparison
scenario_name_2="loving-dingo-vocalize-Horseradish"

## Specify Satellite Data for 1st scenario
sat_data_1="equilibrium"

## Specify x for 2nd scenario
sat_data_2="equilibrium"

## Description of 1st scenario
description_1="Satellite feedback behavior"

## Description of 2nd Scenario
description_2="Open-access behavior"

## Name of save folder
compare_folder="comparison-benchmark-eqm-satfeedback"


### 3. PAIRWISE SCENARIO COMPARISONS

Rscript "compare-two-scenarios.R" $scenario_name_1 $scenario_name_2 $model_type $model_type $sat_data_1 $sat_data_2 $model_horizon $model_horizon "$description_1" "$description_2" $compare_folder

## Tell them you're done buddy!
echo "R scenario comparisons complete!"
# Good job! 

####################
