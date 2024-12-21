import numpy as np
import scipy 
import meep as mp
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src')
from Materials import * 
from FDTD_2Dstructures import geo_2DSlab1L_RHole

import os 

### Resolution 
resolution = 20 

### Boundary layers 
# PML 
dboundary = 2.0  # PML/Absorber thickness 
pml_layers = [mp.PML(direction = mp.X,
                     thickness = dboundary),
              mp.PML(direction = mp.Z,
                     thickness = dboundary)
            ]
 
abs_layers = [mp.Absorber(direction = mp.X,
                          thickness = dboundary),
              mp.Absorber(direction = mp.Z,
                          thickness = dboundary)
            ]

### The number of unit cells along the line y = 0
Ncell = 8

### Padding block
pad = 3.0 

### Geometrical parameters 
d = np.sqrt(2.0) # The diagonal of one square unit cell
h = 0.35   # Thickness of the layer 
b = 0.41   # The edge length of the undeformed square hole
e = 0.1    # The deformation parameter 

# The total size of the structure along the x-axis
structurex = Ncell*d 

# The total size of the structure along the y-axis
structurey = d 

# The height of the unit cell along the z-direction (not including the PML/Absorber layer)
Lz = 3*h 

### Size of the simulation cell 
sx = structurex + 2*(pad + dboundary)
sy = d      # no PML, no pad 
sz = Lz + 2*dboundary 

### Define the simulation cell 
cell = mp.Vector3(sx,sy,sz)

### The materials 
Mater = Si 
Envir = PMMA 

### The vertices of the unit cells and rhombus holes
vertice_cell = [
    mp.Vector3(0.5*d,0,0),
    mp.Vector3(0,0.5*d,0),
    mp.Vector3(-0.5*d,0,0),
    mp.Vector3(0,-0.5*d,0)
]

vertice = [
    mp.Vector3(b*(1+e)/(1-e)/np.sqrt(2),0,0),
    mp.Vector3(0,b*(1-e)/(1+e)/np.sqrt(2),0),
    mp.Vector3(-b*(1+e)/(1-e)/np.sqrt(2),0,0),
    mp.Vector3(0,-b*(1-e)/(1+e)/np.sqrt(2),0)
]

### The source 
fcen = 0.255    # pulse center frequency 
df = 0.015      # pulse width 
nfreq = 501     # number of frequencies
component = mp.Ey   # the component 
sources = [
    mp.Source(
        mp.ContinuousSource(frequency=fcen,is_integrated=True),
        component = component,
        center = mp.Vector3(-0.5*sx+dboundary+0.5*pad,0,0),
        size = mp.Vector3(0,structurey,h)
    )
]

# The array of frequencies
freq_array = np.linspace(fcen-df,fcen+df,nfreq)

##### ===================================================================================
##### Compute the flux spectrum
freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0),
                     size = mp.Vector3(0,structurey,h),
                     direction = mp.X)

##### ====================================================================================
##### The position of the monitor
pt = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0)

##### =================================================================================
##### GEOMETRY 
geometry = geo_2DSlab1L_RHole(d,h,b,e,vertice,Mater,Envir,Ncell,sx,structurey)

##### ==================================================================================
##### Define the simulation 
sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = abs_layers,
        geometry = geometry,
        sources = sources,
        k_point = mp.Vector3(0,0,0),
	    ensure_periodicity = True,
        resolution = resolution
    )

##### ==================================================================================
##### Transmitted flux
trans = sim.add_flux(fcen,df,nfreq,freg)

##### ==================================================================================
##### Run the simulation
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ey", mp.at_every(0.5, mp.output_efield_y)),
        until=200)