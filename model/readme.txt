# README FILE FOR LANDLAB MODELING
#
# PURPOSE: Contains python scripts needed to simulate aspect dependent soil production
#
# INPUT FILES:
# 'swath_hillside_main.py' is the main Landlab script
# 'swath_hillside_functions.py' contains functions used in main script
#
# SAVED FILES:
# 'out_initial_state.csv' is initial x-position [m], surface z [m], bedrock z [m] of nodes
# 'out_final_state.csv' is final x-position [m], surface z [m], bedrock z [m] of nodes
# 'out_time_series.csv' is time [yrs], mean elevation [m], mean soil depth [m] of nodes
# 'output_images/***.png' are save maps from the model run
#
# DESCRIPTION:
# This model script includes process descriptions for depth-dependent, nonlinear creep and 
# humped soil production. The model domain is a quasi-1D swath with open boundaries on the
# left and right boundaries. Random perturbations to topography and soil depth are set as 
# a fraction of the surface by the user. Soil transport is simulated using the component
# 'DepthDependentTaylorDiffuser' in Landlab. Soil production is simulated using a modified 
# version of the soil production function proposed by Strudley et al. (2008) that includes
# a slope correction to represent surface-normal weathering. Saved .csv file can be uses to 
# probe model outputs. Maps with transects are also saved as .png files in a directory 
# called 'output images' for making animated .gifs.
#
# REFERENCES:
# Strudley, M. W., Murray, A. B., & Haff, P. K. (2006). Emergence of pediments, tors, and 
# piedmont junctions from a bedrock weathering–regolith thickness feedback. Geology, 34(10),
# 805-808.
