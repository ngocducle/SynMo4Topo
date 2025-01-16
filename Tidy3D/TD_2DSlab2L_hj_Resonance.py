import numpy as np
import scipy 
import xarray as xr 
import matplotlib.pyplot as plt

import tidy3d as td
from tidy3d import web 
from tidy3d.plugins.resonance import ResonanceFinder

import sys
sys.path.insert(0,'../src')
from Tidy3D_2Dstructures import td_2DSlab2L_RHoleP_hj

# Geometrical parameters 
a = 1               # Period (micrometer)
d = a*np.sqrt(2)    # Diagonal length
h = 0.35*a          # Slab thickness 
b1 = 0.35*a         # Edge of the undeformed square hole (left-hand side)
e1 = 0.05           # deformation parameter (left-hand side)
b2 = 0.41*a         # Edge of the undeformed square hole (right-hand side)
e2 = 0.05           # deformation parameter (right-hand side)
dist = 0.1          # Distance between 2 slabs
delta = 0.50        # Relative displacement (fraction of d)

# Number of unit cells for each side along the diagonal direction (set to be x)
Ncell = 28

# Wavelength and frequency range of interest 
# C_0 = speed of light (micrometer / s)
freq_range = np.array((0.24,0.27)) # unit of c/a
freq_scale = (td.constants.C_0/a)  # frequency scale determined by the lattice constant (Hz)
freq_range_Hz = (freq_range[0]*freq_scale,freq_range[1]*freq_scale) # frequency range (Hz)
lambda_array = (a/freq_range[1],a/freq_range[0])
print(lambda_array)

# Gaussian pulse parameters
freq0 = np.sum(freq_range_Hz)/2             # central frequency (Hz)
freqw = 0.3*(freq_range_Hz[1]-freq_range_Hz[0])  # pulse width (Hz) 
lambda0 = td.constants.C_0 / freq0       # central wavelength (micrometer)
print(f"Central wavelength lambda0 = {lambda0} micrometer")

# Runtime
runtime_fwidth = 1000.0  # in units of 1/frequency bandwidth of the source
t_start_fwidth = 5.0    # time to start monitoring after source has decayed, units of 1/frequency bandwidth
run_time = runtime_fwidth / freq0 
print(f"Total runtime = {(run_time*1e12):.2f} ps")
t_start = t_start_fwidth / freq0 

# Frequencies and wavelengths of monitor
Nfreq = 1001 
monitor_freqs = np.linspace(freq_range[0],freq_range[1],Nfreq) # unit: c/a
monitor_lambda = td.constants.C_0 / monitor_freqs # unit: micrometer

# Materials 
ep_slab = 3.54**2   # Si
mat_slab = td.Medium(permittivity=ep_slab,name='silicon')
ep_envir = 1.46**2  # PMMA/SiO2
mat_envir = td.Medium(permittivity=ep_envir,name='PMMA')

# Space between PhC slabs and PML
pad_dimless = 1.5
pad = pad_dimless * lambda_array[-1] # unit: micrometer

# Grid size (micrometer)
dl = lambda_array[0]/30/np.sqrt(ep_slab)   # 30 grids per smallest wavelength in medium (unit: micrometer)
print(f"dl = {dl*1000} nm")

# Simulation size 
sim_size = Lx,Ly,Lz = (2*pad+2*Ncell*d,d,2*h+dist+2*lambda_array[-1])

# Simulation structure
sim_structures = td_2DSlab2L_RHoleP_hj(d,h,b1,e1,b2,e2,Lx,pad,dist,delta,mat_envir,mat_slab,Ncell)

### Boundary conditions
bspec = td.BoundarySpec(
    x = td.Boundary(minus=td.PML(num_layers=12),plus=td.PML(num_layers=12)),
    y = td.Boundary(minus=td.Periodic(),plus=td.Periodic()),
    z = td.Boundary(minus=td.PML(num_layers=12),plus=td.PML(num_layers=12)),
)

### Point dipole source, gaussian pulse, placed at the interface
"""source = td.PointDipole(
    source_time = td.GaussianPulse(
        freq0 = freq0,
        fwidth = freqw,
    ),
    center = (0,0,0),
    polarization = 'Ey',
    name = "dipole",
)"""

rng = np.random.default_rng(12345)

num_dipoles = 3

dipole_positions = rng.uniform(
    [0,-0.5*d,-0.5*dist-0.5*h],[0,0.5*d,0.5*dist+0.5*h],[num_dipoles,3]
)

dipole_phases = rng.uniform(0,2*np.pi,num_dipoles)

pulses = []
dipoles = []
for i in range(num_dipoles):
    pulse = td.GaussianPulse(freq0=freq0,fwidth=freqw,phase=dipole_phases[i])
    pulses.append(pulse)
    dipoles.append(
            td.PointDipole(
                source_time = pulse,
                center = tuple(dipole_positions[i]),
                polarization = 'Ey',
                name = "dipole_"+str(i),
                )
            )

### Field time monitor
"""monitor = td.FieldTimeMonitor(
    fields = ["Ey"],
    center = (0,0,0.5*dist+0.3*h),
    size = (0,0,0),
    start = t_start,
    name = "monitor_time",
)"""

num_monitors = 5

monitor_positions = rng.uniform(
    [0,-0.5*d,-0.5*dist-0.5*h],[0,0.5*d,0.5*dist+0.5*h],[num_monitors,3]
)

monitors = []
for i in range(num_monitors):
    monitors.append(
        td.FieldTimeMonitor(
            fields = ["Ey"],
            center = tuple(monitor_positions[i]),
            size = (0,0,0),
            start = t_start,
            name = "monitor_time_"+str(i),
    )
)


##### Define the simulation 
sim = td.Simulation(
    center = (0,0,0),
    size = sim_size,
    grid_spec = td.GridSpec.uniform(dl=dl),
    structures = sim_structures,
    sources = dipoles,
    monitors = monitors,
    run_time = run_time,
    shutoff = 0,
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

fig,ax = plt.subplots(1,2,tight_layout=True,figsize=(10,6))
sim.plot(x=0.0,ax=ax[0])
sim.plot(y=0.0,ax=ax[1])
plt.savefig(namesave+'Structure_yz_xz_planes.png')
plt.show()

# Initializes job, puts task on server (but doesnt run it)
job = web.Job(simulation=sim, task_name="job", verbose=True)

# Estimate the maximum cost
estimated_cost = web.estimate_cost(job.task_id)

##### Running simulation
sim_data = web.run(sim,task_name='2DSlab2L_transmission',path=f'data/2DSlab2L-Ncell_{Ncell:d}-delta_{delta:.4f}.hdf5')