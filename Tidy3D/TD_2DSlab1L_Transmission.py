import numpy as np
import tidy3d as td 
from tidy3d import web
from tidy3d.constants import C_0
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src')
from Tidy3D_2Dstructures import td_2DSlab1L_RHoleP

# Geometrical parameters  
a = 0.78         # Period (micrometer)
d = a*np.sqrt(2) # Diagonal length 
h = 0.35*a       # Slab thickness
b = 0.35*a       # Edge of the undeformed square hole
e = 0.1          # deformation parameter

# Wavelength and frequency range
freq_range_dimless = (0.24,0.27) # Dimensionless unit fa/c
freq_range = (freq_range_dimless[0]*C_0/a,freq_range_dimless[1]*C_0/a)
#print(freq_range)
lambda_range = (a/freq_range_dimless[1],a/freq_range_dimless[0])
freq0 = np.sum(freq_range)/2

# Central frequency, frequency pulse width 
lambda0 = C_0 / freq0 
width = 0.3 
freqw = width*(freq_range[1]-freq_range[0])

# Runtime 
runtime = 100.0 # in unit of 1/frequency bandwidth of the source 
t_stop = runtime/freqw 
print(f"Total runtim <= {t_stop*1e12} ps")

# Number of unit cell along the diagonal direction (set to be x)
Ncell = 7 

# Materials
ep_slab = 3.54**2  # Si 
mat_slab = td.Medium(permittivity=ep_slab,name='silicon')
ep_envir = 1.46**2 # PMMA/SiO2
mat_envir = td.Medium(permittivity=ep_envir,name='PMMA')

# Space between PhC slabs and PML
pad_dimless = 1.0 
pad = pad_dimless * lambda_range[-1]

# Grid size # micrometer 
dl = lambda_range[0]/30/np.sqrt(ep_slab) # 30 grids per smallest wavelength in medium
print(f"dl = {dl*1000} nm")

# Simulation size 
sim_size = Lx,Ly,Lz = (2*pad+Ncell*d,d,h+2*h*lambda_range[-1])

"""### The environment
envir = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,td.inf),
            ),
        medium = mat_envir,
        name = 'envir',
        )

### Define the slab
slab = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,h),
            ),
        medium = mat_slab,
        name = 'slab',
        )

### Initially, the simulated structure contains the environment and the slab
sim_structures = [envir,slab]"""

### Define the rhombus hole 
"""vertices = np.array([ (b*(1+e)/(1-e)/np.sqrt(2),  0),
                    (0,  b*(1-e)/(1+e)/np.sqrt(2)),
                    (-b*(1+e)/(1-e)/np.sqrt(2), 0),
                    (0, -b*(1-e)/(1+e)/np.sqrt(2)) ])"""

"""hole = td.Structure(
        geometry = td.PolySlab(
            axis = 2,
            reference_plane = 'middle',
            slab_bounds = [-0.5*h,0.5*h],
            vertices = vertices,
            ),
        medium = mat_envir,
        name = 'hole',
        )"""

#sim_structures.append(hole)

"""Nhalf = int((Ncell-1)/2)
print('Nhalf = '+str(Nhalf))

# Central line
for j in np.arange(-Nhalf,Nhalf+1):
    holename = 'holec_'+str(j)

    vertices = np.array(
            [
                (j*d + b*(1+e)/(1-e)/np.sqrt(2),0),
                (j*d, b*(1-e)/(1+e)/np.sqrt(2)),
                (j*d - b*(1+e)/(1-e)/np.sqrt(2),0),
                (j*d, -b*(1-e)/(1+e)/np.sqrt(2))
                ]
            )

    hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

    sim_structures.append(hole)

# Upper line
for j in np.arange(-Nhalf,Nhalf):
    holename = 'holeu_'+str(j)

    vertices = np.array(
            [
                ( (j+0.5)*d + b*(1+e)/(1-e)/np.sqrt(2),0.5*d),
                ( (j+0.5)*d, b*(1-e)/(1+e)/np.sqrt(2) +0.5*d),
                ( (j+0.5)*d - b*(1+e)/(1-e)/np.sqrt(2),0.5*d),
                ( (j+0.5)*d, -b*(1-e)/(1+e)/np.sqrt(2)+0.5*d)
                ]
            )

    hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

    sim_structures.append(hole)

# Lower line 
for j in np.arange(-Nhalf,Nhalf):
    holename = 'holel_'+str(j)

    vertices = np.array(
            [
                ( (j+0.5)*d + b*(1+e)/(1-e)/np.sqrt(2),-0.5*d),
                ( (j+0.5)*d, b*(1-e)/(1+e)/np.sqrt(2) -0.5*d),
                ( (j+0.5)*d - b*(1+e)/(1-e)/np.sqrt(2),-0.5*d),
                ( (j+0.5)*d, -b*(1-e)/(1+e)/np.sqrt(2)-0.5*d)
                ]
            )

    hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

    sim_structures.append(hole)"""

sim_structures = td_2DSlab1L_RHoleP(d,b,h,e,mat_envir,mat_slab,Ncell)

### Boundary conditions
bspec = td.BoundarySpec(
        x = td.Boundary(minus=td.PML(num_layers=12),plus=td.PML(num_layers=12)),
        y = td.Boundary(minus=td.Periodic(),plus=td.Periodic()),
        z = td.Boundary(minus=td.PML(num_layers=12),plus=td.PML(num_layers=12)),
        )

##### Define the simulation
sim = td.Simulation(
        center = (0,0,0),
        size = sim_size,
        grid_spec = td.GridSpec.uniform(dl=dl),
        structures = sim_structures,
        sources = [],
        monitors = [],
        run_time = t_stop,
        shutoff = 1e-7,
        boundary_spec = bspec,
        normalize_index = None,
        )


##### Plot the structure
fig,ax = plt.subplots(1,2,tight_layout=True,figsize=(10,4))
sim.plot(z=0.0,ax=ax[0])
sim.plot(x=0.0,ax=ax[1])
plt.show()
