import numpy as np
import tidy3d as td 
from tidy3d.constants import C_0 
import matplotlib.pyplot as plt 

### Geometrical parameters 
a = 1   # Period (micrometer)
freq_range_dimless = (0.24,0.27)    # Dimensionless unit fa/c

### Import data from file
sim_data_norm = td.SimulationData.from_file('data/2DSlab2L-normalization.hdf5')
sim_data_050 = td.SimulationData.from_file('data/2DSlab2L-Ncell_8-delta_0.5000.hdf5')
sim_data_051 = td.SimulationData.from_file('data/2DSlab2L-Ncell_8-delta_0.5100.hdf5')
sim_data_052 = td.SimulationData.from_file('data/2DSlab2L-Ncell_8-delta_0.5200.hdf5')
sim_data_053 = td.SimulationData.from_file('data/2DSlab2L-Ncell_8-delta_0.5300.hdf5')
sim_data_054 = td.SimulationData.from_file('data/2DSlab2L-Ncell_8-delta_0.5400.hdf5')

### Retrieve the power flux through the monitor plane 
transmission_norm = sim_data_norm['flux'].flux

transmission050 = sim_data_050['flux'].flux 
transmission050 = np.array(transmission050 / transmission_norm)
transmission051 = sim_data_051['flux'].flux 
transmission051 = np.array(transmission051 / transmission_norm)
transmission052 = sim_data_052['flux'].flux 
transmission052 = np.array(transmission052 / transmission_norm)
transmission053 = sim_data_053['flux'].flux 
transmission053 = np.array(transmission053 / transmission_norm)
transmission054 = sim_data_054['flux'].flux 
transmission054 = np.array(transmission051 / transmission_norm)

print(transmission052)

Trans = np.array([transmission050,transmission051,transmission052,transmission053])
print(Trans)

##### Plot the figure
fig,ax = plt.subplots(1,figsize=(10,8),tight_layout=True)
plt.plot(transmission050,label='delta = 0.50')
plt.plot(transmission051,label='delta = 0.51')
plt.plot(transmission052,label='delta = 0.52')
plt.plot(transmission053,label='delta = 0.53')
#plt.plot(transmission054)
ax.set_xlim(400,600)
plt.legend()
plt.show()

##### Plot the spectrum
fig,ax = plt.subplots()
im = ax.pcolormesh(Trans.T,shading='gouraud',cmap='hot')
fig.colorbar(im,ax=ax)
plt.show()