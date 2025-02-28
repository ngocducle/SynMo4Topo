import numpy as np
import tidy3d as td
from tidy3d import web 
from tidy3d.constants import C_0
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src')
from Tidy3D_2Dstructures import td_2DSlab2L_RHoleP_hj 

# Geometrical parameters 
a = 1             # Period (micrometer)
d = a*np.sqrt(2)  # Diagonal length 
h = 0.35*a        # Slab thickness 
b1 = 0.46*a       # Edge of the undeformed square hole (left-hand side)
e1 = 0.1          # deformation parameter (left-hand side)
b2 = 0.30*a       # Edge of the undeformed square hole (right-hand side)
e2 = 0.1          # deformation parameter (right-hand side)
dist = 0.1        # Distance between 2 slabs
delta = 0.50      # Relative displacement (fraction of d)

# Number of unit cells for each side along the diagonal direction (set to be x)
Ncell = 8

# Wavelength and frequency range 
freq_range_dimless = (0.245,0.265) # Dimensionless unit fa/c
freq_range = (freq_range_dimless[0]*C_0/a,freq_range_dimless[1]*C_0/a)
lambda_range = (a/freq_range_dimless[1],a/freq_range_dimless[0])
freq0 = np.sum(freq_range)/2

# Central frequency, frequency pulse width 
lambda0 = C_0 / freq0 
width = 0.3 
freqw = width*(freq_range[1]-freq_range[0])
print(f"Central wavelength = {1000*lambda0} nm")

# Runtime 
runtime = 50.0 # in unit of 1/frequency bandwidth of the source
t_stop = runtime / freqw 
print(f"Total runtime <= {t_stop*1e12} ps")

# Frequencies and wavelengths of monitor 
Nfreq = 1001 
monitor_freqs = np.linspace(freq_range[0],freq_range[1],Nfreq)
monitor_lambda = C_0 / monitor_freqs

# Materials 
ep_slab = 3.54**2 # Si 
mat_slab = td.Medium(permittivity=ep_slab,name='silicon')
ep_envir = 1.46**2 # PMMA/SiO2
mat_envir = td.Medium(permittivity=ep_envir,name='PMMA')

# Space between PhC slabs and PML 
pad_dimless = 1.5
pad = pad_dimless * lambda_range[-1]

# Grid size (micrometer)
dl = lambda_range[0]/30/np.sqrt(ep_slab) # 30 grids per smallest wavelength in medium
print(f"dl = {dl*1000} nm")

# Simulation size 
sim_size = Lx,Ly,Lz = (2*pad+2*Ncell*d,d,2*h+dist+2*lambda_range[-1])

# Simulation structure 
sim_structures = td_2DSlab2L_RHoleP_hj(d,h,b1,e1,b2,e2,Lx,pad,dist,delta,mat_envir,mat_slab,Ncell)

### Boundary conditions
bspec = td.BoundarySpec(
    x = td.Boundary(minus=td.PML(num_layers=12),plus=td.PML(num_layers=12)),
    y = td.Boundary(minus=td.Periodic(),plus=td.Periodic()),
    z = td.Boundary(minus=td.PML(num_layers=12),plus=td.PML(num_layers=12)),
)

##### Plane wave source, gaussian pulse, placed just in advance of the slab (towards positive x)
source = td.PlaneWave(
    source_time = td.GaussianPulse(
        freq0 = freq0,
        fwidth = freqw,
    ),
    size = (0,td.inf,2*h+dist),
    center = (-0.5*Lx+0.5*pad,0,0),
    direction = '+',
    angle_theta = 0,
    pol_angle = 0,
    name = 'planewave'
)

##### Flux monitor 
monitor = td.FluxMonitor(
    center = (0.5*Lx-0.5*pad,0,0),
    size = (0,td.inf,2*h+dist),
    freqs = monitor_freqs,
    name = 'flux',
)

##### Time monitor 
monitor_time = td.FieldTimeMonitor(
    center = (0.5*Lx-0.5*pad,0,0),
    size = (0,td.inf,h),
    interval = 1000,
    fields = ['Ex','Ey','Ez','Hx','Hy','Hz'],
    name = 'time',
)

##### Define the simulation 
sim = td.Simulation(
    center = (0,0,0),
    size = sim_size,
    grid_spec = td.GridSpec.uniform(dl=dl),
    structures = sim_structures,
    sources = [source],
    monitors = [monitor],
    run_time = t_stop,
    shutoff = 1e-5,
    boundary_spec = bspec,
    normalize_index = None,
)

##### Plot the structure 
namesave = f'Ncell_{Ncell:d}-'

fig,ax = plt.subplots(3,1,tight_layout=True,figsize=(10,6))
sim.plot(z=0.5*(h+dist),ax=ax[0])
sim.plot(z=0.0,ax=ax[1])
sim.plot(z=-0.5*(h+dist),ax=ax[2])
plt.savefig(namesave+'Structure_xy_plane.png')
plt.show()

fig,ax = plt.subplots(1,1,tight_layout=True,figsize=(10,8))
sim.plot(x=0.0)
plt.savefig(namesave+'Structure_yz_planes.png')
plt.show()

fig,ax = plt.subplots(1,1,tight_layout=True,figsize=(20,8))
sim.plot(y=0.0)
plt.savefig(namesave+'Structure_xz_planes.png')
plt.show()

##### Check probe and source
f,(ax1,ax2) = plt.subplots(1,2,tight_layout=True,figsize=(8,4))
plot_time = 2e-12
ax1 = sim.sources[0].source_time.plot(times=np.linspace(0,plot_time,1001),ax=ax1)
ax1.set_xlim(0,plot_time)
ax1.legend(('source amplitude',))
ax2 = sim.sources[0].source_time.plot_spectrum(times=np.linspace(0,sim.run_time,10001),val='abs',ax=ax2)
fill_max = 6e-14
ymax = 6e-14
ax2.fill_between(freq_range,[-0e-16,-0e-16],[fill_max,fill_max],alpha=0.4,color='g')
ax2.legend(('source spectrum','measurement'))
ax2.set_ylim(-1e-16,ymax)
plt.savefig(namesave+'Sources_profile.png')
plt.show()

##### Running simulation
#sim_data = web.run(sim,task_name='2DSlab2L_transmission',path=f'data/2DSlab2L-Ncell_{Ncell:d}-delta_{delta:.4f}.hdf5')