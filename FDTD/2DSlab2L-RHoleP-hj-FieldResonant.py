import numpy as np
import scipy 
import meep as mp 
import matplotlib.pyplot as plt 

import sys
sys.path.insert(0,'../src')
from Materials import * 
from FDTD_2Dstructures import geo_2DSlab2L_RHoleP_hj_sameMater

import os 

### Resolution 
resolution = 16

### Boundary layers
dboundary = 3.0  # PML/Absorber thickness 
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

mix_layers = [mp.Absorber(direction = mp.X,
                          thickness = dboundary),
              mp.PML(direction = mp.Z,
                    thickness = dboundary)
            ]

### The source
fcen = 0.255    # pulse center frequency 
df = 0.015      # pulse width 
nfreq = 501     # number of frequencies

# The array of frequencies and wavelength
freq_array = np.linspace(fcen-df,fcen+df,nfreq)
wvl_array = 1/freq_array 

### The value of the synthetic momentum
q = 0.0

### The number of unit cells along the line y = 0
Ncell = 8

### Padding block
pad = 2.0*wvl_array[0] 

### Geometrical parameters 
d = np.sqrt(2.0)   # The diagonal of one square unit cell
h = 0.35           # Thickness of the two layers

# The layer 1 
b1 = 0.41   # The edge length of the undeformed square hole
e1 = 0.1    # The deformation parameter 

# The layer 2 
b2 = 0.35   # The edge length of the undeformed square hole 
e2 = 0.1    # The deformation parameter 

# The interlayer distance 
dist = 0.1

# The total size of the structure along the x-axis
structurex = 2*Ncell*d 

# The toal size of the structure along the y-axis
structurey = d 

# The total thickness of the bilayer (z-direction)
hbilayer = 2*h + dist 

# The height of the unit cell along the z-direction (not including the PML/Absorber layer)
Lz = 3*hbilayer  

### Size of the simulation cell 
sx = structurex + 2*(pad + dboundary)
sy = d      # no PML, no pad, PBC 
sz = Lz + 2*dboundary

### Convert to shift
delta = d*(q + 0.5) 

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

vertice1 = [
    mp.Vector3(b1*(1+e1)/(1-e1)/np.sqrt(2),0,0),
    mp.Vector3(0,b1*(1-e1)/(1+e1)/np.sqrt(2),0),
    mp.Vector3(-b1*(1+e1)/(1-e1)/np.sqrt(2),0,0),
    mp.Vector3(0,-b1*(1-e1)/(1+e1)/np.sqrt(2),0)
]

vertice2 = [
    mp.Vector3(b2*(1+e2)/(1-e2)/np.sqrt(2),0,0),
    mp.Vector3(0,b2*(1-e2)/(1+e2)/np.sqrt(2),0),
    mp.Vector3(-b2*(1+e2)/(1-e2)/np.sqrt(2),0,0),
    mp.Vector3(0,-b2*(1-e2)/(1+e2)/np.sqrt(2),0)
]

##### The source
component = mp.Ey   # the component 
sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df,is_integrated=True),
        component = component,
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(0,structurey,hbilayer)
    )
]

##### ==================================================================================
##### GEOMETRY
geometry = geo_2DSlab2L_RHoleP_hj_sameMater(d,h,hbilayer,delta,
                                            vertice1,vertice2,
                                            Mater,Envir,Ncell,sx,sy)

##### ====================================================================================
##### The position of the monitor
pt = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0)

##### ==================================================================================
##### Define the simulation 
sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = abs_layers,
        geometry = geometry,
        sources = sources,
        k_point = mp.Vector3(0,0,0), # PBC
	    ensure_periodicity = True,   # PBC
        resolution = resolution
    )

##### ==================================================================================
##### Run the simulation 
sim.run(
        # Output the dielectric profile at the beginning
        mp.at_beginning(mp.output_epsilon),
        # Output the field Ey at every 1 time unit
        mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0,0,0), fcen, df)),
        # Finish the simulation when the field at the point pt decays,
        # check the maximum of the field every dt
        until_after_sources = 400
        )

sim.run(mp.at_every(1/fcen/20, mp.output_efield_y), until=1/fcen)