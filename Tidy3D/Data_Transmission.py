import numpy as np
import tidy3d as td 
import matplotlib.pyplot as plt

### Import data from file
sim_data = td.SimulationData.from_file('data/data_2DSlab1L_Ncell6_pad4.875.hdf5')
sim_data0 = td.SimulationData.from_file('data/data0_2DSlab1L_Ncell0_pad4.875.hdf5')

### Retrieve the power flux through the monitor plane 
transmission0 = sim_data0['flux'].flux
transmission = sim_data['flux'].flux
transmission_normalized = transmission / transmission0 

### Plot the figure
fig,(ax1,ax2) = plt.subplots(2,1,figsize=(4,7),tight_layout=True)
transmission0.plot(ax=ax1,label='no slab')
#transmission.plot(ax=ax1,label='with slab')
transmission_normalized.plot(ax=ax2)
ax1.legend()
ax1.set_title('raw transmission')
ax2.set_title('normalized transmission')
plt.show()
