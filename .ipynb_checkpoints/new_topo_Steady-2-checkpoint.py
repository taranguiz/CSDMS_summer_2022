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

# Model dimensions
xmax= 3000
ymax= 1000
dxy= 10
nrows = int(ymax/dxy)
ncols = int(xmax/dxy)
print(nrows)
print(ncols)

#Instantiate model grid
mg= RasterModelGrid ((nrows,ncols), dxy)
#add field topographic elevation
mg.add_zeros("node", "topographic__elevation")

np.random.seed(seed=5000)
#creating initial model topography
random_noise = (np.random.rand(len(mg.node_y)))

#add the topo to the field
mg["node"]["topographic__elevation"] +=random_noise

# add field 'soil__depth' to the grid
mg.add_zeros("node", "soil__depth", clobber=True)

# Set  5m of initial soil depth at core nodes
mg.at_node["soil__depth"][mg.core_nodes] = 5.0  # meters

# Add field 'bedrock__elevation' to the grid
mg.add_zeros("bedrock__elevation", at="node")
# Sum 'soil__depth' and 'bedrock__elevation'
# to yield 'topographic elevation'
mg.at_node["bedrock__elevation"][:] = mg.at_node["topographic__elevation"]
mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]
mg.add_zeros("node", "soil_production__rate", clobber=True)
# soil_production_rate= mg.at_node["soil_production__rate"]

rock= mg.at_node["bedrock__elevation"]
soil=mg.at_node["soil__depth"]
z=mg.at_node["topographic__elevation"]

# Geomorphic parameters
# uplift
uplift_rate= 5 *1e-4

#Hillsope Geomorphology for DDTD component
H=5 # original soil depth
Hstar= 0.1 # characteristic transport depth, m
V0= 0.1 #transport velocity coefficient
D= 0.001 #V0 *Hstar  #effective(maximum) diffusivity

#Fluvial Erosion for SPACE Large Scale Eroder
K_sed=10*1e-5 #sediment erodibility
K_br= 5*1e-5 #bedrock erodibility
F_f=0.5 #fraction of fine sediment
phi= 0.5 #sediment porosity
H_star=Hstar #sediment entrainment lenght scale
Vs= 1 #velocity of sediment
m_sp= 0.5 #exponent ondrainage area stream power
n_sp= 1 #exponent on channel slope in the stream power framework
sp_crit_sed=0 #sediment erosion threshold
sp_crit_br=0 #bedrock erosion threshold


# instantiate components
#Weathering
expweath=ExponentialWeatherer(mg, soil_production__maximum_rate=0.001, soil_production__decay_depth=Hstar)
# Hillslope with Diffuser
ddd=DepthDependentDiffuser(mg, linear_diffusivity=D,
                                  soil_transport_decay_depth=Hstar)

#Flow Router
fr=PriorityFloodFlowRouter(mg, flow_metric='D8', suppress_out=True, runoff_rate=0.5)
#SPACE Large Scale
space= SpaceLargeScaleEroder(mg,
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

mg.set_closed_boundaries_at_grid_edges(
    bottom_is_closed=False,
    left_is_closed=True,
    right_is_closed=True,
    top_is_closed=True,
)

# Now the for loop to do landscape evolution
fig = plt.figure(figsize=[8, 8])
imshow_grid(mg, z, cmap='terrain', grid_units=['m', 'm'])
# plt.show()

#timing
tmax=2000000
dt=100
model_time=np.arange(0,tmax,dt)
iterations=len(model_time)

for i in range(iterations):
    print(i)
    # z_beg=z[:]
    # expweath.run_one_step()
    expweath.calc_soil_prod_rate()
    ddd.run_one_step(dt)
    fr.run_one_step()
    space.run_one_step(dt)
    
    z[mg.core_nodes] += uplift_rate * dt
    rock[mg.core_nodes] += uplift_rate * dt
    
    z_end=z[:]
    
    print(sum(z_end))
    
    if i%200 ==0:
        fig = plt.figure(figsize=[8, 8])
        imshow_grid(mg, z, cmap='terrain', grid_units=['m', 'm'])
        # mg.save('output_new_topo_ddd_2/topo_%s_yrs.asc' %(int(i * dt)))
        plt.title('Topography after ' + str(int((i * dt))) + ' years')
        plt.show()
        plt.savefig('output_new_topo_ddd_2/topo_%s_yrs.png' % (int(i * dt)), dpi=300, facecolor='white')
mg.save('output_new_topo_ddd_2/finaltopo.asc', at='node')    
    # if float(diff) < float(uplift_rate):
    #     print('steady state reached')
    #     break


