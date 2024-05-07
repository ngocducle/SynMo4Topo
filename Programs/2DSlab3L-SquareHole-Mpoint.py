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
resolution = 32

### Geometrical parameters
h = 0.30      # Thickness of each slab
edge = 0.9    # Edge of the square hole
Lz = 10.0     # The height of the unit cell along the z-direction 
dist = 0.05   # Interlayer distance

#delta_x = 0.0 # Shift along the x-axis
#delta_y = 0.0 # Shift along the y-axis

Ndelta = 51
delta_array = np.linspace(0.0,0.5,Ndelta)

print('# Thickness of the slab h = '+str(h))
print('# Edge of the square hole b = '+str(edge))
print('# The height of the simulation cell Lz = '+str(Lz)) 

### Number of bands 
num_bands = 20 
print('# The number of bands to simulate: '+str(num_bands))

### The title and the name of the files 
namesave = '2DSlab3L-SquareHole-h_'+str(h)+'-b_'+str(edge)

### The array of frequencies of bands
bands = np.zeros((Ndelta,num_bands))
print(np.shape(bands))

###### We run over the array delta_array 
for i in range(Ndelta):

    delta_x = delta_array[i]
    delta_y = 0.0 

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
        mp.Vector3(0.5,0.5,0.0)
    ]

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
    print(np.shape(freqs))

    ### Save the frequencies to the array bands
    bands[i,:] = freqs 

##### Print the arrays delta_array and bands to a file
with open(namesave+'-frequencies.txt','w') as file:
    for i in range(Ndelta):
        file.write('%.8f    ' % delta_array[i])

        for j in range(num_bands):
            file.write('%.8f    ' % bands[i,j]) 

        file.write('\n')





