# dose-finding-simulation
Assess the impact of the choice of dose-finding method in phase 1 on later phase success

data_generating_functions.R: contains functions to generate individual patient data for different types of outcomes (DLT, response, survival)

simulate_phases_functions.R: core functions that users directly interact with. Perform simulation given input specification

combine_results.R: function to combine simulation results from different phases and calculate summary statistics throughout the phases.

main_simulation.qmd: file to demonstrate how to use the functions and present an example of Imatinib.
