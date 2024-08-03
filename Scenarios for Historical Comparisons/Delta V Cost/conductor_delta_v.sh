#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=04:00:00
#SBATCH --partition=amilan

module purge
module load matlab

####################
### 1. INITIAL DEFINITIONS
# Don't need to change, defined here for reducing typo errors later

# Solver script name
MATLAB_SCRIPT="iam_solver.m" 

# Behavior types
launch_pattern_type_equilibrium="equilibrium"
launch_pattern_type_feedback="sat_feedback"

####################
### 2. SCENARIO DEFINITIONS
# Change to set scenarios and parallelization

# Choose propagator
model_type="MOCAT" # Either "MOCAT" or "GMPHD"

# Define the scenario to run. Store them in an array. Should be valid names of parameter set CSV files. 
## See examples in scenarios/parsets and compare to files named --parameters.csv for how to create new ones.
scenario_files=(
                "scenarios/parsets/delta_v_0.csv"
                "scenarios/parsets/delta_v_2000.csv"
                )

# Define the length of the horizon over which to solve (years)
model_horizon=35

# Set number of workers for parallelization
n_workers=40

####################
### 3.  EXECUTION
# For each scenario and launch behavior: 1. Run the MATLAB script; 2. Run the R script, loop over scenarios with the scenario_file inputs; 3. Run comparisons.

# Step 3.1: Run simulations

# Initialize empty array to store unique names for step 3.2
unique_name_array=()

# Main simulation loop
for scenario_file in "${scenario_files[@]}"
do
    # Call the scenarioNamer function with the parameter file as an argument
    matlab -nosplash -nodesktop -r "scenarioNamer('$scenario_file', '$model_type'); exit;"

    # Read the unique name from the scenario_name.txt file
    unique_name=$(cat scenario_name.txt)
    unique_name_array+=("$unique_name")

    # Launch behavior 1: Equilibrium, MOCAT or GMPHD
    matlab -nosplash -nodesktop -r "iam_solver('$unique_name', '$model_type', '$launch_pattern_type_equilibrium', '$scenario_file', '$model_horizon', '$n_workers'); exit;"

    # Launch behavior 2: Feedback, MOCAT only
    matlab -nosplash -nodesktop -r "iam_solver('$unique_name', 'MOCAT', '$launch_pattern_type_feedback', '$scenario_file', '$model_horizon', '$n_workers'); exit;"

    # Tell them you're done buddy!
    echo "Scenario $scenario_file complete!"
done