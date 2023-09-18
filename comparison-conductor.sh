#!/bin/bash

### 1. INITIAL DEFINITIONS

# Define R command line arguments for scenario parameters
## Change as needed
model_type="MOCAT" # Either "MOCAT" or "GMPHD"
## Don't need to change, defined here for reducing typo errors later
launch_pattern_type_equilibrium="equilibrium"
launch_pattern_type_feedback="sat_feedback"

# Define the R script(s) to run
R_SCRIPT="analytics.R" 

### 2. SCENARIO DEFINITIONS
# Select the scenarios to process. Store them in an array. Directly use human-readable names here, this script only runs R files.
scenario_names=(
                "polluting-flamingo-enunciate-Cilantro"
                "philosophical-arctic-fox-bends-LemonBalm"
                "measuring-goose-add-Knapweed"
                )

# Define the length of the horizon over which to solve (years)
model_horizon=35

### 3.  EXECUTION
# Loop over scenarios with the scenario_file inputs and run the analytics script. Can skip this if they've already been run in conductor.sh.
for scenario_name in "${scenario_names[@]}"
do
    # Launch behavior 1: Equilibrium, MOCAT or GMPHD
    Rscript "analytics.R" $scenario_name $model_type $launch_pattern_type_equilibrium $model_horizon

    # Launch behavior 2: Feedback, MOCAT only
    Rscript "analytics.R" $scenario_name $model_type $launch_pattern_type_feedback $model_horizon

    # Tell them you're done buddy!
    echo "Scenario ~~ $scenario_name ~~ processing complete!"
done

# Run comparisons. This piece will not have been run in conductor.sh.

# ## Benchmark: equilibrium vs sat feedback
mkdir scenarios/comparison--benchmark-eqm-satfeedback

Rscript "compare-two-scenarios.R" "${scenario_names[1]}" "${scenario_names[1]}" "MOCAT" "MOCAT" "sat_feedback" "equilibrium" "35" "35" "Satellite feedback behavior" "Open-access behavior" "comparison--benchmark-eqm-satfeedback"

# Rscript "compare-two-scenarios.R" "polluting-flamingo-enunciate-Cilantro" "polluting-flamingo-enunciate-Cilantro" "MOCAT" "MOCAT" "sat_feedback" "equilibrium" "35" "35" "Satellite feedback behavior" "Open-access behavior" "comparison--benchmark-eqm-satfeedback"

## 5 vs 25 year eol comparison
mkdir scenarios/comparison--5-vs-25-yr-eol

Rscript "compare-two-scenarios.R" "${scenario_names[2]}" "${scenario_names[1]}" "MOCAT" "MOCAT" "equilibrium" "equilibrium" "35" "35" "25-year disposal rule" "5-year disposal rule" "comparison--5-vs-25-yr-eol"

## 5 year no-ouf vs 25 year ouf comparison
mkdir scenarios/comparison--5-year-no-ouf-vs-25-year-ouf

Rscript "compare-two-scenarios.R" "${scenario_names[1]}" "${scenario_names[3]}" "MOCAT" "MOCAT" "equilibrium" "equilibrium" "35" "35" "5-year disposal, w/o OUF" "25-year disposal, w/ OUF" "comparison--5-year-no-ouf-vs-25-year-ouf"

## Tell them you're done buddy!
echo "Scenario comparisons complete!"
# Good job! 
