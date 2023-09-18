#!/bin/bash

### 1. INITIAL DEFINITIONS

# Define the MATLAB script to run
MATLAB_SCRIPT="iam_solver.m" # Don't need to change, defined here for reducing typo errors later

# Define MATLAB command line arguments for scenario parameters
## Change as needed
model_type="MOCAT" # Either "MOCAT" or "GMPHD"
## Don't need to change, defined here for reducing typo errors later
launch_pattern_type_equilibrium="equilibrium"
launch_pattern_type_feedback="sat_feedback"

### 2. SCENARIO DEFINITIONS
# Define the scenario to run. Store them in an array.
scenario_files=(
                "benchmark"
                "scenarios/parsets/25yr-rule-parset.csv"
                "scenarios/parsets/25yr-rule-tax-0.5-parset.csv"
                )

# Initialize empty array to store unique names
unique_name_array=()

# Define the length of the horizon over which to solve (years)
model_horizon=10

# Set number of workers for parallelization
n_workers=40

### 3.  EXECUTION
# For each scenario and launch behavior: 1. Run the MATLAB script; 2. Run the R script. Loop over scenarios with the scenario_file inputs. 3. Run comparisons.
# Step 1
for scenario_file in "${scenario_files[@]}"
do
    # Call the scenarioNamer function with the parameter file as an argument
    matlab -nodisplay -nosplash -nodesktop -r "scenarioNamer('$scenario_file', '$model_type'); exit;"

    # Read the unique name from the scenario_name.txt file
    unique_name=$(cat scenario_name.txt)
    unique_name_array+=("$unique_name")

    # Launch behavior 1: Equilibrium, MOCAT or GMPHD
    matlab -nodisplay -nosplash -nodesktop -r "iam_solver('$unique_name', '$model_type', '$launch_pattern_type_equilibrium', '$scenario_file', '$model_horizon', '$n_workers'); exit;"
    # Rscript "analytics.R"  $unique_name $model_type $launch_pattern_type_equilibrium $model_horizon

    # Launch behavior 2: Feedback, MOCAT only
    matlab -nodisplay -nosplash -nodesktop -r "iam_solver('$unique_name', 'MOCAT', '$launch_pattern_type_feedback', '$scenario_file', '$model_horizon', '$n_workers'); exit;"
    Rscript "analytics.R"  $unique_name $model_type $launch_pattern_type_feedback $model_horizon

    # Tell them you're done buddy!
    echo "Scenario $scenario_file complete!"
done

# Step 2
for scenario_name in "${unique_name_array[@]}"
do
    # Launch behavior 1: Equilibrium, MOCAT or GMPHD
    Rscript "analytics.R" $scenario_name $model_type $launch_pattern_type_equilibrium $model_horizon

    # Launch behavior 2: Feedback, MOCAT only
    Rscript "analytics.R" $scenario_name $model_type $launch_pattern_type_feedback $model_horizon

    # Tell them you're done buddy!
    echo "Scenario ~~ $scenario_name ~~ processing complete!"
done

## Tell them you're done buddy!
echo "R script execution complete!"
# Good job! 

# Step 3

# ## Benchmark: equilibrium vs sat feedback
mkdir scenarios/comparison--benchmark-eqm-satfeedback

Rscript "compare-two-scenarios.R" "${unique_name_array[1]}" "${unique_name_array[1]}" "MOCAT" "MOCAT" "sat_feedback" "equilibrium" "35" "35" "Satellite feedback behavior" "Open-access behavior" "comparison--benchmark-eqm-satfeedback"

# Rscript "compare-two-scenarios.R" "polluting-flamingo-enunciate-Cilantro" "polluting-flamingo-enunciate-Cilantro" "MOCAT" "MOCAT" "sat_feedback" "equilibrium" "35" "35" "Satellite feedback behavior" "Open-access behavior" "comparison--benchmark-eqm-satfeedback"

## 5 vs 25 year eol comparison
mkdir scenarios/comparison--5-vs-25-yr-eol

Rscript "compare-two-scenarios.R" "${unique_name_array[2]}" "${unique_name_array[1]}" "MOCAT" "MOCAT" "equilibrium" "equilibrium" "35" "35" "25-year disposal rule" "5-year disposal rule" "comparison--5-vs-25-yr-eol"

## 5 year no-ouf vs 25 year ouf comparison
mkdir scenarios/comparison--5-year-no-ouf-vs-25-year-ouf

Rscript "compare-two-scenarios.R" "${unique_name_array[1]}" "${unique_name_array[3]}" "MOCAT" "MOCAT" "equilibrium" "equilibrium" "35" "35" "5-year disposal, w/o OUF" "25-year disposal, w/ OUF" "comparison--5-year-no-ouf-vs-25-year-ouf"

## Tell them you're done buddy!
echo "Scenario comparisons complete!"
# Good job! 


