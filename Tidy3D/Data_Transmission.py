import numpy as np
import tidy3d as td 
from tidy3d.constants import C_0
import matplotlib.pyplot as plt

### Geometrical parameters 
a = 1       # Period (micrometer)
freq_range_dimless = (0.24,0.27) # Dimensionless unit fa/c

### Import data from file
sim_data = td.SimulationData.from_file('data/2DSlab2L-Ncell_8-delta_0.5000.hdf5')
sim_data0 = td.SimulationData.from_file('data/2DSlab2L-normalization.hdf5')

### Retrieve the power flux through the monitor plane 
transmission0 = sim_data0['flux'].flux
transmission = sim_data['flux'].flux
transmission_normalized = transmission / transmission0 

### Plot the figure
fig,(ax1,ax2) = plt.subplots(2,1,figsize=(4,7),tight_layout=True)
transmission0.plot(ax=ax1,label='no slab')
transmission.plot(ax=ax1,label='with slab')
transmission_normalized.plot(ax=ax2)
ax1.legend()
ax1.set_title('raw transmission')
ax2.set_title('normalized transmission')
plt.savefig('transmission-Ncell_8.png')
plt.show()