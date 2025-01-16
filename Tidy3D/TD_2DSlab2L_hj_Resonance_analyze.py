import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import tidy3d as td 
from tidy3d import web
from tidy3d.plugins.resonance import ResonanceFinder

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
Ncell = 8
namesave = f'Ncell_{Ncell:d}-'

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
runtime_fwidth = 200.0  # in units of 1/frequency bandwidth of the source
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

### Import data from file
sim_data = td.SimulationData.from_file('data-time_1000/2DSlab2L-Ncell_8-delta_0.5000.hdf5')

### Plot the field with respect to time
fig = plt.figure()
plt.plot(sim_data.monitor_data["monitor_time_4"].Ey.t,
         np.real(sim_data.monitor_data["monitor_time_4"].Ey.squeeze()))
plt.title("FieldTimeMonitor data")
plt.xlabel("t")
plt.ylabel("Ey")
plt.show()

### Fourier transform the data to see the resonances at the band frequencies
field = sim_data.monitor_data["monitor_time_4"].Ey.squeeze().real

fmesh = np.fft.fftfreq(field.size, np.mean(np.diff(field.t)))
spectrum = np.fft.fft(field)

mask = (fmesh > freq_range_Hz[0]) & (fmesh < freq_range_Hz[1])

plt.plot(
        fmesh[mask],
        np.abs(spectrum)[mask],
        )
plt.title("Spectrum at single wavevector")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.savefig(namesave+'spectrum.png')
plt.show()

# Use the ResonanceFinder plugin to find the band frequencies
resonance_finder = ResonanceFinder(freq_window=tuple(freq_range))
resonance_data = resonance_finder.run(signals=sim_data.data)
modes = resonance_data.to_dataframe()
print(modes)