# -*- coding: utf-8 -*-
"""
Functions used by 'swath_hillside_main.py'

Purpose:
    Make Landlab grids. 
    Modify Lanlab grids.
    Plot results.


Requirements:
    numpy, landlab, matplotlib, seaborn
"""

## PACKAGES ##
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from landlab import RasterModelGrid as RMG

## DEFINE FUNCTIONS ##

## INITIALIZE LANDLAB GRID ##
# Make initial swath profile with open boundaries along sides
def make_initial_grid(xy ,dx, e_rate, K_creep, sp0, h_star):
    grid = RMG(xy, xy_spacing=dx)                           # make Landlab grid
    centerx = (np.round(xy[1]/2))*dx                        # swath center
    zmax = (e_rate / (2 * K_creep)) * (((xy[1]*dx)/2)**2)   # steady state relief for linear creep [m]
    ssh = -1 * (h_star * np.log(e_rate / sp0))              # steady state soil depth for exp wxing [m]
    
    print(f'Initial Soil Depth: {ssh:0.2f} m')

    # ADD topography field
    z = grid.add_zeros('topographic__elevation', at='node') # make elevation
    z[:] = zmax - (e_rate/(2*K_creep))*((grid.x_of_node-centerx)**2)     # steady state topogrpahy

    # ADD soil depth field
    sd = grid.add_zeros('soil__depth', at='node')          # make soil depth field
    sd[:] = ssh                                           # initial soil depth (background)

    # DERIVE bedrock elevation   
    zb = grid.add_zeros('bedrock__elevation', at='node')    # make bedrock z
    zb[:] = z[:] - sd[:]                                    # set bedrock z

    # ADD empty fields
    sp = grid.add_zeros("soil_production__rate", at="node")  # soil production rate

    # SET boundary conditions
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)

    return grid, z, zb, sd, sp

## CALCULATION ON LANDLAB GRIDS ##
#  Calculate soil production using coside correction
def soil_production_geometry(grid, sp, sd, SP0, d1, d2, k1):
    sl = grid.calc_slope_at_node(elevs='topographic__elevation', method='patch_mean', ignore_closed_nodes=True, return_components=False)
    sp_hum = SP0 * (np.exp((-1)*sd[grid.core_nodes]/d1)-(k1*np.exp((-1)*sd[grid.core_nodes]/d2)))*(1/np.cos(sl[grid.core_nodes]))# humped fun
    sp[grid.core_nodes] = sp_hum
    
# Modify creep coefficient as a function of local aspect
def modify_diffusion(grid,z,cr,K_creep,asp):    
    [low, high] = [(4/5)*K_creep, (6/5)*K_creep]  # high a low values (2X contrast)
    [a, b, c, d] = [45, 135, 225, 310]            # define bounds for classes (degrees)
        
    cr[(asp>a) & (asp<b)] = low                   # low values facing 'right' (N-facing)
    cr[(asp>c) & (asp<d)] = high                  # high values facing 'left' (S-facing)

# Modify soil production rates as a function of local aspect
def modify_production(grid,sp,asp):
    [a, b, c, d] = [45, 135, 225, 310]            # define bounds for classes (degrees)
        
    sp[(asp>a) & (asp<b)] *= (6/5)              # high values values facing 'right' (N-facing)
    sp[(asp>c) & (asp<d)] *= (4/5)               # low values facing 'left' (S-facing)
            
## PLOTTING ##
# Plot current state of Landlab grid
def plot_state(z,zb,sd,grid,xy,dx,fname):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, height_ratios=[1,1,5]) 
    
    IBM = ['#648FFF','#785EF0','#DC267F','#FE6100','#FFB000']   # IBM colorblind friendly
    
    # MAP of surface elevation
    draw1 = np.reshape(z, xy)
    #ax1 = sns.heatmap(draw1, vmin=0, vmax=np.max(z[:]), xticklabels=False, 
    #                  yticklabels=False, ax=ax1, cmap='viridis', cbar=False, square=True)
    ax1 = sns.heatmap(draw1, vmin=0, vmax=150, xticklabels=False, 
                      yticklabels=False, ax=ax1, cmap='cividis', cbar=False, square=True)
    ax1.set_title(f'Mean Elevation: {np.mean(z[grid.core_nodes]):0.1f} m', fontdict=dict(fontsize=24))
    
    # MAP of soil depths
    draw2 = np.reshape(sd, xy) 
    #ax2 = sns.heatmap(draw2, vmin=0, vmax=np.max(sd[:]), xticklabels=False, 
    #                  yticklabels=False, ax=ax2, cmap='viridis', cbar=False, square=True)
    ax2 = sns.heatmap(draw2, vmin=0, vmax=2, xticklabels=False, 
                      yticklabels=False, ax=ax2, cmap='cividis', cbar=False, square=True)
    ax2.set_title(f'Mean Soil Depth: {np.mean(sd[grid.core_nodes]):0.2f} m', fontdict=dict(fontsize=24))

    # PLOT of swath profile
    x=grid.x_of_node
    
    
    # central five transects
    cen = xy[1]*int(xy[0]/2)      # rounds down
    add = xy[1]
    
    ax3.plot(x[cen-(2*add):cen-add-1], z[cen-(2*add):cen-add-1], c=IBM[2], lw=1.5, label=None)  # 2 above
    ax3.plot(x[cen-add:cen-1], z[cen-add:cen-1], c=IBM[2], lw=1.5, label=None)  # above
    ax3.plot(x[cen:cen+add-1], z[cen:cen+add-1], c=IBM[2], lw=1.5, label=None)  # central
    ax3.plot(x[cen+add:cen+(2*add)-1], z[cen+add:cen+(2*add)-1], c=IBM[2], lw=1.5, label=None)  # below
    ax3.plot(x[cen+(2*add):cen+(3*add)-1], z[cen+(2*add):cen+(3*add)-1], c=IBM[2], lw=1.5, label=None)  # 2 below

    ax3.set(xlim=(0, (xy[1]-1)*dx), ylim=(-1,175))
    ax3.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                         bottom=True, top=True, left=True, right=True, direction='in')
    ax3.tick_params(axis='both', which='major', labelsize=16)
    
    ax3.set_xlabel('Distance [m]', size=24) 
    ax3.set_ylabel('Elevation [m]', size=24)
    #ax3.legend(loc='upper right', frameon=False, fontsize=24)
    
    fig.set_size_inches(8,9)
    plt.savefig(fname, bbox_inches='tight')
    plt.show()

# Plot time evolution soil depth and relief
def plot_evolution(t,z,sd,z_ss,sd_ss):
    fig, ax1 = plt.subplots(1,1)

    color = 'r' 
    ax1.plot(t, z, color)
    ax1.plot([0,1.1*np.max(t)], [z_ss, z_ss], c=color, ls=':')
    ax1.set(xlim=(0, np.max(t)), ylim=(0,1.1*np.max(z)))
    ax1.set_xlabel('Time') 
    ax1.set_ylabel('Relief (m)', color=color, labelpad=4)
    ax1.xaxis.set_tick_params(labelcolor='k', pad=10, length=0)
    ax1.yaxis.set_tick_params(labelcolor=color, pad=10, length=0)

    ax2 = ax1.twinx()  # twin axis
    color = 'b'
    ax2.plot(t, sd, color=color)
    ax2.plot([0,1.1*np.max(t)], [sd_ss, sd_ss], c=color, ls=':')
    ax2.set(xlim=(0, np.max(t)), ylim=(0,1.1*np.max(sd)))
    ax2.set_ylabel('Mean Soil Depth (m)', color=color, rotation=270, labelpad=25) 
    ax2.yaxis.set_tick_params(labelcolor=color, pad=10, length=0)
    
    fig.set_size_inches(8,3)
    plt.savefig('output2.png', bbox_inches='tight')
    plt.show()
    