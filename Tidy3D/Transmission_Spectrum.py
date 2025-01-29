import numpy as np
import tidy3d as td 
from tidy3d.constants import C_0 
import matplotlib.pyplot as plt 

### Geometrical parameters 
a = 1   # Period (micrometer)
freq_range_dimless = (0.24,0.26)    # Dimensionless unit fa/c

### Import data from file
sim_data_norm = td.SimulationData.from_file('data-0.240-0.260/2DSlab2L-normalization.hdf5')
sim_data_040 = td.SimulationData.from_file('data-0.240-0.260/2DSlab2L-Ncell_8-delta_0.4000.hdf5')

### Retrieve the power flux through the monitor plane 
transmission_norm = sim_data_norm['flux'].flux

transmission040 = sim_data_040['flux'].flux 
transmission040 = np.array(transmission040 / transmission_norm)

#Trans = np.array([transmission050,transmission051,transmission052,transmission053])
#print(Trans)

##### Plot the figure
fig,ax = plt.subplots(1,figsize=(10,8),tight_layout=True)
plt.plot(transmission040,label='delta = 0.40')
#plt.plot(transmission051,label='delta = 0.51')
#plt.plot(transmission052,label='delta = 0.52')
#plt.plot(transmission053,label='delta = 0.53')
#plt.plot(transmission054)
#ax.set_xlim(400,600)
plt.legend()
plt.savefig('Transmission_0.240_0.260.png')
plt.show()

##### Plot the spectrum
#fig,ax = plt.subplots()
#im = ax.pcolormesh(Trans.T,shading='gouraud',cmap='hot')
#fig.colorbar(im,ax=ax)
#plt.show()