
#!/usr/bin/env python3
# import time
import numpy as np
import matplotlib.pyplot as plt

#from Landlab
from landlab import RasterModelGrid, imshow_grid, imshowhs_grid
from landlab.io import read_esri_ascii

#Hillslope geomorphology
from landlab.components import ExponentialWeatherer
from landlab.components import DepthDependentTaylorDiffuser
from landlab.components import DepthDependentDiffuser

#Fluvial Geomorphology and Flow routing
from landlab.components import FlowDirectorMFD #trying the FlowDirectorMFD
from landlab.components import FlowAccumulator, Space, FastscapeEroder, PriorityFloodFlowRouter
from landlab.components.space import SpaceLargeScaleEroder

# Reading the original topo
(grid, z) = read_esri_ascii('topo142x109.asc',name='topographic__elevation')
grid.set_closed_boundaries_at_grid_edges(bottom_is_closed=False, left_is_closed=True, right_is_closed=True, top_is_closed=True)

xmax = 142 # x grid dimension in meters
ymax = 109 # y grid dimension in meters
dxy = 1.0 # grid step in meters
nrows = int(ymax/dxy)
ncols = int(xmax/dxy)

# Choosing the fault location
fault_loc = 120
fault_nodes = np.where(grid.node_y==fault_loc)[0]

# plotting parameters (grid plotting & initial conditions)
figsize = [8,8] # size of grid plots
shrink = 0.35 # amount of colorbar shrinkage for plots (0-1). 1 = not shrunk. 0 = nonexistent.
limits = [0,20] # elevation limits for grid plots

fig = plt.figure(figsize=figsize)
imshow_grid(grid,z,grid_units=['m','m'], cmap='gray', shrink=shrink)
plt.title('Original Topography')
plt.savefig('output_ddtd/Original Topography',dpi=300,facecolor='white')
#plt.show()

# uplift
uplift_rate= 1 *1e-4

#Hillsope Geomorphology for DDTD component
H=10 # original soil depth
Sc= 0.7 #critical slope
Hstar= 0.1 # characteristic transport depth, m
V0= 0.01 #transport velocity coefficient
D= V0 *Hstar  #effective(maximum) diffusivity

#Fluvial Erosion for SPACE Large Scale Eroder
K_sed=0.001 #sediment erodibility
K_br= 0.00001 #bedrock erodibility
F_f=0.5 #fraction of fine sediment
phi= 0.5 #sediment porosity
H_star=Hstar #sediment entrainment lenght scale
Vs= 0.01 #velocity of sediment
m_sp= 0.5 #exponent ondrainage area stream power
n_sp= 1 #exponent on channel slope in the stream power framework
sp_crit_sed=0 #sediment erosion threshold
sp_crit_br=0 #bedrock erosion threshold


#Creating fields that components require
grid.add_zeros("node", "soil__depth", clobber=True) #add field to the grid
grid.at_node["soil__depth"]=grid.at_node["soil__depth"]+H
grid.at_node["bedrock__elevation"]=grid.at_node["topographic__elevation"] - grid.at_node["soil__depth"]
grid.add_zeros("node", "soil_production__rate", clobber=True)
soil_production_rate= grid.at_node["soil_production__rate"]
rock= grid.at_node["bedrock__elevation"]
soil=grid.at_node["soil__depth"]
# figsize = [16,4] # size of grid plots
# fig, ax = plt.subplots(figsize=figsize)
# x = grid.node_x[fault_nodes]
soil_level = rock + soil
# ax.plot(x, soil_level[fault_nodes], 'orange', linewidth=2, markersize=12, label='soil')
# ax.plot(x, rock[fault_nodes], linewidth=2, markersize=12, label='bedrock')
# plt.title('Original cross-Profile topography at fault location')
# ax.set_xlabel('X (m)')
# ax.set_ylabel('Depth (m)')
# ax.legend(loc='lower right')
#plt.show()

#timing
tmax=100000
dt=10
model_time=np.arange(0,tmax,dt)
iterations=len(model_time)

# instantiate components
# Hillslope with Taylor Diffuser
ddtd=DepthDependentTaylorDiffuser(grid,
                                  soil_transport_velocity=V0,
                                  soil_transport_decay_depth=Hstar,
                                  nterms=2,
                                  dynamic_dt=True,
                                  if_unstable='warn')
#Flow Router
fr=PriorityFloodFlowRouter(grid, flow_metric='D8', suppress_out=True)
#SPACE Large Scale
space= SpaceLargeScaleEroder(grid,
                             K_sed=K_sed,
                             K_br=K_br,
                            F_f=F_f,
                            phi=phi,
                            H_star=Hstar,
                            v_s=Vs,
                            m_sp=m_sp,
                            n_sp=n_sp,
                            sp_crit_sed=0,
                             sp_crit_br=0)
# Now the for loop to do landscape evolution

for i in range(iterations):
    print('before landscape processes')
    z[grid.core_nodes] += uplift_rate *dt
    rock[grid.core_nodes] += uplift_rate * dt
    ddtd.run_one_step(dt)
    #ddd.run_one_step(dt)
    fr.run_one_step()
    space.run_one_step(dt)
    #erosion_rate=SpaceLargeScaleEroder._calc_erosion_rates()
    if i%100  == 0:
        fig = plt.figure(figsize=[8, 8])
        imshow_grid(grid, z, cmap='gray', grid_units=['m', 'm'], shrink=shrink)
        plt.title('Topography after '+str(int((i*dt)))+' years')
        plt.savefig('output_ddtd/topo_%s_yrs.png' % (int((i * dt)+dt)), dpi=300, facecolor='white')
    print(i)
    #print(erosion_rate)
    print('after landscape processes')
#plt.show()
print('job done')
#
# fig = plt.figure(figsize=[8, 8])
# imshow_grid(grid, z, cmap='gray', grid_units=['m', 'm'], shrink=shrink)
# plt.title('Topography after '+str(int(i*dt +dt))+' years')
# plt.show()
