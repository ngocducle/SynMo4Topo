import numpy as np 
import scipy 
import meep as mp 
from meep import mpb 

import sys
sys.path.insert(0,'../src/')
from Materials import * 
from DielectricProfile import * 
from ExportData import * 
from BandStructure import * 

### Resolution 
resolution = 8

### Geometrical parameters
h = 0.30      # Thickness of each slab
edge = 0.9    # Edge of the square hole
Lz = 10.0     # The height of the unit cell along the z-direction 
dist = 0.05   # Interlayer distance
delta_x = 0.0 # Shift along the x-axis
delta_y = 0.0 # Shift along the y-axis

print('# Thickness of the slab h = '+str(h))
print('# Edge of the square hole b = '+str(edge))
print('# The height of the simulation cell Lz = '+str(Lz)) 

### Number of bands 
num_bands = 20 
print('# The number of bands to simulate: '+str(num_bands))

### Number of points to interpolate the k-points
Nk = 19

### Define the lattice
geometry_lattice = mp.Lattice(
    size = mp.Vector3(1.0,1.0,Lz),
    basis1 = mp.Vector3(1.0,0.0),
    basis2 = mp.Vector3(0.0,1.0)
)

### Define the materials 
Mater = Si 
Envir = PMMA 

### Define the geometry
geometry = [
    mp.Block(
        center = mp.Vector3(0.0,0.0,0.0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir 
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,0.0),
        size = mp.Vector3(1.0,1.0,1.8*h),
        material = Mater
    ),

    mp.Block(
        center = mp.Vector3(delta_x,delta_y,0.0),
        size = mp.Vector3(edge,edge,1.8*h),
        material = Envir 
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,1.4*h),
        size = mp.Vector3(1.0,1.0,h),
        material = Mater 
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,1.4*h),
        size = mp.Vector3(0.3*edge,0.3*edge,h),
        material = Envir 
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,-1.4*h),
        size = mp.Vector3(1.0,1.0,h),
        material = Mater 
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,-1.4*h),
        size = mp.Vector3(0.3*edge,0.3*edge,h),
        material = Envir 
    )
]

### The k_points 
k_points = [
    mp.Vector3(0.0,0.0,0.0),    # Gamma 
    mp.Vector3(0.5,0.0,0.0),    # X 
    mp.Vector3(0.5,0.5,0.0),    # M 
    mp.Vector3(0.0,0.0,0.0)     # Gamma 
]

k_points = mp.interpolate(Nk,k_points)

### The ModeSolver
ms = mpb.ModeSolver(
    geometry = geometry,
    geometry_lattice = geometry_lattice,
    k_points = k_points,
    resolution = resolution,
    num_bands = num_bands
)

### Run the simulation 
ms.run_zeven()

### Extract the frequencies of the modes from the ModeSolver
freqs = ms.all_freqs

### The number of elements in k_points
#number = np.arange(len(ms.k_points))

### The title and the name of the files 
namesave = '2DSlab3L-SquareHole-h_'+str(h)+'-b_'+str(edge)

### Output the dielectric profile with Ncellx x Ncelly unit cells 
### for zmin <= z <= zmax, with Nz values of z 
Ncellx = 5
Ncelly = 5
zmin = -2.0*h
zmax = 2.0*h 
Nx = 300
Ny = 300
Nz = 21

# Calculate the dielectric profile in planes parallel to Oxy
x_plot,y_plot,z_array,epsilon_xy_array \
    = DielectricProfileXY(ms,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz)

# Plot the dielectric profile, the name of the figure is:
#           namesave+'-z_'+str(k)+'.png'
# where k is the number of the value of z in the array z_array
PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave,'Yes')

### Print the band structure to file
PrintBandStructure(freqs,number,namesave)

PlotBand_BrillouinZone(number,freqs,Nk,namesave,'Yes')