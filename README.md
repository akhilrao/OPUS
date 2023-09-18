# OPUS (Orbital Debris Propagators Unified with Economic Systems): An Integrated Assessment Model for Satellites and Orbital Debris

This repository provides a comprehensive package for the replication and extension of the Integrated Assessment Model (IAM) described in the accompanying paper. The IAM focuses on simulating the evolution of orbital environments and is implemented primarily using MATLAB for the core simulation and R for post-simulation analytics.

## Repository Structure

### Main Scripts and Directories

- **`conductor.sh`**: This Bash shell script serves as the main driver for running IAM. It orchestrates the setup, execution, and post-processing of multiple scenarios.

- **`iam_solver.m`**: This is the main solver file implemented in MATLAB. It performs the time-stepped simulation of the orbital environment based on a set of input parameters and scenarios.

### Supporting Scripts

- **`analytics.R`**: This R script handles the analysis and visualization of the simulation results.

- **`MOCAT4S/` and `GMPHD/`**: These directories contain the MATLAB scripts for the MOCAT and GMPHD propagators respectively.

- **`historical-data/`**: This folder contains the historical data needed for various scenarios.

- **`scenarios/`**: This folder contains CSV files that define different scenarios for the model.

### Additional Resources

- **`scenarioNamer.m`**: Generates unique, human-readable names for different scenarios based on the model parameters.

- **`modifyParameters.m` and `set_econ_parameters.m`**: These MATLAB scripts help in parameter modification and initialization respectively.

## Workflow

### Initialization and Scenario Setup (`conductor.sh`)

1. **Initialization**: Sets up initial parameters like MATLAB solver to use (`iam_solver.m`), optional analytics script (`analytics.R`), and the type of propagator (MOCAT or GMPHD).

2. **Scenario Setup**: Initializes scenarios to run and sets up paths for associated CSV files.

3. **Execution and Post-Processing Loop**: Executes the IAM for each scenario and runs the analytics.

### IAM Solver (`iam_solver.m`)

This MATLAB script includes several key steps, such as initialization, parameter modification, model propagation, and result saving. For an in-depth understanding, consult Algorithm 1 in the accompanying paper.

## Usage

1. Clone the repository.
2. Navigate to the root directory.
3. Run `./conductor.sh` in the shell.

## Acknowledgments

This work was supported by NASA ROSES.
