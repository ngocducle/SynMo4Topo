import numpy as np
import scipy 
import meep as mp 
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src')
from Materials import * 
from FDTD_2Dstructures import geo_2DSlab2L_RHoleP_hj_sameMater

import os 

##### ===============================================================================
##### Calculate the Ey field profile through a heterojunction 
##### knowing the frequency of the electric field mode 
##### The output data are analyzed by the file 2DSlab2L-PlotFields.py
##### ===============================================================================

### Resolution 
resolution = 20

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

### The number of unit cells along the line y = 0
Ncell = 8

### Padding block
pad = 3 

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

# The interlayer shift
delta = 0.5

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
        mp.GaussianSource(fcen,fwidth=df),
        component = component,
        center = mp.Vector3(-delta/2-0.027,+0.026,0.0),
        size = mp.Vector3(0,0,0)
    )
]

 ##### GEOMETRY
geometry = geo_2DSlab2L_RHoleP_hj_sameMater(d,h,hbilayer,delta,
                                            vertice1,vertice2,
                                            Mater,Envir,Ncell,sx,sy)

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

##### Volume to extract the field
Field_vol = mp.Volume(
    center = mp.Vector3(0,0,0),
    size = mp.Vector3(structurex,structurey,Lz)
)

##### The array of x,y,z-direction
x_array = np.linspace(-Ncell*d,Ncell*d,int(structurex*resolution)+2)
y_array = np.linspace(-0.5*d,0.5*d,int(d*resolution)+2)
z_array = np.linspace(-0.5*Lz,0.5*Lz,int(Lz*resolution)+2)

##### The discrete Fourier transform
dft_obj = sim.add_dft_fields([mp.Ey],fcen,0,1,where=Field_vol)


##### Run the simulation
sim.run(
    until_after_sources=mp.stop_when_fields_decayed(
        300,
        mp.Ey,
        mp.Vector3(-d-0.2,0.18),
        1e-7)
    )

##### Get the dielectric and electric field profiles
eps_data = sim.get_array(vol=Field_vol,component=mp.Dielectric)
Ey_data = sim.get_dft_array(dft_obj,mp.Ey,0)

##### Save the x_array, y_array, z_array 
print('Shape of x_array:')
print(np.shape(x_array))
print('Shape of y_array:')
print(np.shape(y_array))
print('Shape of z_array:')
print(np.shape(z_array))
print('Shape of eps_data:')
print(np.shape(eps_data))
print('Shape of Ey_data:')
print(np.shape(Ey_data))

np.savetxt('x_array.txt',x_array)
np.savetxt('y_array.txt',y_array)
np.savetxt('z_array.txt',z_array)

### Save the dielectric profile along z
for iz in np.arange(len(z_array)): 
    namesave_eps = f'dielectric-z_{iz:d}.txt'
    np.savetxt(namesave_eps,eps_data[:,:,iz])
    
    namesave_Ey = f'Field_Ey-z_{iz:d}.txt'
    np.savetxt(namesave_Ey,Ey_data[:,:,iz])

### Save the integrated electric field Ey 
np.savetxt('dft_field_along_x.txt',
           np.column_stack(
               (
                   x_array,
                   np.sum(
                       np.square(np.absolute(Ey_data)),
                       axis=(1,2) # take integration over dydz
                    ) 
                )
            )
        )