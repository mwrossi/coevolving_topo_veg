# -*- coding: utf-8 -*-
"""
Script using Landlab modeling library to simulate hillside swath profile evolution.

Purpose:
    Produce patchy soils with tors. 
    Simulate asymmetric soil production rates.

Two main process laws:
    (1) Soil creep uses DepthDependentTaylorDiffuser
    (2) Soil production uses a hybrid of exponential and humped functions.

Requirements:
    numpy, time, landlab, math, matplotlib, seaborn, scipy
    swath_hillside_functions.py    
"""

## IMPORT PACKAGES AND LIBRARIES ##
import os
import numpy as np
import time as tic
import swath_hillside_functions as fun
from landlab.components import DepthDependentTaylorDiffuser

## SETUP PARAMETERS ##
# Model Domain
xy = [21,200]      # grid size; shoot for 1:10 for nice plotting
dx = 2             # node spacing [m]
dt = 1             # time step [yrs]
tt = 20000        # total time [yrs]
nt = int(tt//dt)   # int number time steps

# Process Parameters
U = 0.000075       # uplift rate [m/yr]
K = 0.01           # creep coefficient [m2/yr]
h_star = 0.5       # creep depth scaling [m]
s_crit = 0.80      # creep critical slope
SP0 = 0.000240      # production max for exponential [m/yr]
d1 = 0.5           # production depth scaling [m]; 0.5 in Strudley et al. (2006)
d2 = 0.1           # production depth scaling [m]; 0.03 in Strudley et al. (2006)
k1 = 0.9          # production ratio at zero [-]; 0.8 in Strudley et al. (2006)

# Calculate steady state values
ss_depth = -1 * (d1 * np.log(U / SP0))           # ss depth for exp wxing [m]
ss_relief = (U / (2 * K)) * (((xy[1]*dx)/2)**2)  # ss relief for lin creep [m]
print(f'Steady State Soil Depth: {ss_depth:0.2f} m')
print(f'Steady State Relief (linear creep): {ss_relief:0.1f} m')

# Make directory to store files
directory = 'output_images'
os.makedirs(directory, exist_ok=True)

## MAKE LANDLAB GRID ##
# Make steady state profile for prior uplift rate
U_factor = 0.5
sd_init = -1 * (d1 * np.log((U*U_factor) / SP0))           # ss depth for exp wxing [m]


# Make parabolic profile with initial soil depth based on prior steady state
[grid,z,zb,sd,sp] = fun.make_initial_grid(xy=xy,dx=dx,e_rate=U*U_factor,K_creep=K,sp0=SP0,h_star=d1)

# Make objects
hgt = 3           # object height [m]
frac = 0.1        # fraction of surface with objects
    
# Place objects on topography
rng = np.random.default_rng(seed=50)
rand = rng.random(z[grid.core_nodes].size)   # random values from 0-1
rand[rand<(1-frac)] = 0                      # no object
rand[rand>=(1-frac)] = 1                     # object
    
z[grid.core_nodes] += rand*hgt

# Calc slope to set initial soil depths
sl = grid.calc_slope_at_node(elevs='topographic__elevation', method='patch_mean', ignore_closed_nodes=False, return_components=False)

# Adjust initial soil depth to range from zero to maximum value
sd[grid.core_nodes] -= sd_init*sl[grid.core_nodes] / np.max(sl[grid.core_nodes])  # scale depth to local gradient

# Recalculate bedrock elevations
zb[grid.core_nodes] = z[grid.core_nodes] - sd[grid.core_nodes]                    # reset bedrock z


## INITIALIZE LANDLAB COMPONENTS ##
# DDTD notes: K = v * hstar; soil depth is modified both by transport and production
diffuse = DepthDependentTaylorDiffuser(grid, slope_crit=s_crit, 
                                       soil_transport_decay_depth=h_star,
                                       dynamic_dt=True, if_unstable='warn',
                                       soil_transport_velocity=K/h_star)

## RUN MODEL ##
st = tic.time()     # start clock

# Lists for storing model states
time = [0]                                  # store time for plots [yrs]
relief = [np.max(z[grid.core_nodes])]       # store relief for plots
soil_depth = [np.mean(sd[grid.core_nodes])] # store mean soil depth for plots
divide = [(xy[1]/2)*dx]                     # store 

# Plot initial condition
fun.plot_state(z,zb,sd,grid,xy,dx,'output_images/000.png')
np.savetxt('out_initial_state.csv', np.c_[grid.x_of_node[grid.core_nodes],zb[grid.core_nodes],z[grid.core_nodes]], delimiter=',')

# LOOP to iterate over processes at time step dt
for i in range(nt):
    fun.soil_production_geometry(grid, sp, sd, SP0, d1, d2, k1)
    asp = grid.calc_aspect_at_node(elevs=z, unit='degrees')    
    fun.modify_production(grid,sp,asp)
    z[grid.core_nodes] += U*dt
    zb[grid.core_nodes] += U*dt 
    diffuse.run_one_step(dt)

    if (i+1) % 1000 == 0:
        time.append(i)
        relief.append(np.max(z[grid.core_nodes]))
        soil_depth.append(np.mean(sd[grid.core_nodes]))
        
        high = grid.x_of_node[z==np.max(z[grid.core_nodes])]
        divide.append(high)
    
    if (i+1) % 10000 == 0:
        if ((i+1)/10000) < 10:
            fname = f'output_images/00{int((i+1)/10000)}.png'
        elif ((i+1)/10000) >= 10 and ((i+1)/10000) < 100:
            fname = f'output_images/0{int((i+1)/10000)}.png'
        else:
            fname = f'output_images/{int((i+1)/10000)}.png'
            
        fun.plot_state(z,zb,sd,grid,xy,dx,fname)
        print(f'Execution time to {(i+1)*dt} years: {(tic.time()-st)/60:0.2f} minutes') 


## SAVE DATA TO .CSV FILE ##
np.savetxt('out_time_series.csv', np.c_[time, relief, soil_depth], delimiter=',')
np.savetxt('out_final_state.csv', np.c_[grid.x_of_node[grid.core_nodes],zb[grid.core_nodes],z[grid.core_nodes]], delimiter=',')
