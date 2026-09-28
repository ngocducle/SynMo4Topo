import numpy as np
import meep as mp
from meep import mpb
import math 
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0,'../src/')
#from DielectricProfile import *
#from BandStructure import * 
#from ExportData import * 

### Light polarization
polarization = 'zeven'
print('# Polarization = '+polarization)

### Resolution 
resolution = mp.Vector3(32,32,32)    # pixels/a
print('# The resolution:'+str(resolution))

### Number of bands 
num_bands = 15
print('# The number of bands to simulate: '+str(num_bands))

### Geometrical parameters 
h = 0.35    # Thickness of the slab
Lz = 5.0    # The height of the unit cell along the z-direction
b = 0.38    # The edge of the square hole

### Number of k-points to interpolate between 2 high-symmetry points
Nk = 29
print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

### Materials 
Mater = mp.Medium(index = 3.54)
Envir = mp.Medium(index = 1.46)

### Define the lattice
alpha = np.pi/4

geometry_lattice = mp.Lattice(
    size = mp.Vector3(1,1,Lz),
    basis_size = mp.Vector3(1,1,1),
    basis1 = mp.Vector3(np.cos(alpha),np.sin(alpha)),
    basis2 = mp.Vector3(np.cos(alpha),-np.sin(alpha))
    #basis1=mp.Vector3(math.sqrt(3)/2, 0.5),
    #basis2=mp.Vector3(math.sqrt(3)/2, -0.5)
)

### Define the geometry
geometry = [
    mp.Block(
        center = mp.Vector3(0.0,0.0,0.0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,0.0),
        size = mp.Vector3(mp.inf,mp.inf,h),
        material = Mater
    ),

    # Note that the vertices of the prism are given in the basis 
    # of primitive lattice vector, so they have the same values 
    # as for the case where the unit cell is not rotated
    mp.Prism(
        vertices = [
            mp.Vector3( b/2, b/2, 0.0 ),
            mp.Vector3( -b/2, b/2, 0.0 ),
            mp.Vector3(-b/2, -b/2, 0.0 ),
            mp.Vector3( b/2, -b/2, 0.0 )
        ],
        height = h,
        axis = mp.Vector3(0.0,0.0,1.0),
        center = mp.Vector3(0.0,0.0,0.0),
        material = Envir 
    )
]

### The k-points
k_points = [
    mp.Vector3(0.0,0.0,0.0),
    mp.Vector3(0.5,0.0,0.0),
    mp.Vector3(0.5,0.5,0.0),
    mp.Vector3(0.0,0.0,0.0)
]

k_points = mp.interpolate(Nk,k_points)

### Mode solver
ms = mpb.ModeSolver(
    geometry = geometry,
    geometry_lattice = geometry_lattice,
    k_points = k_points,
    resolution = resolution,
    num_bands = num_bands
)

### Run the simulation
ms.run_zeven()

### Display the dielectric profile
"""md = mpb.MPBData(rectify=True, periods=3, resolution=50)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)

plt.imshow(converted_eps, interpolation='spline36', cmap='binary')
plt.axis('on')
plt.show()"""



### Extract the frequencies of the modes from the ModeSolver
freqs = ms.all_freqs

### The number of elements in k_points
number = np.arange(len(ms.k_points))

### The title and the name of the files 
namesave = '2DSlab1L-SquareHole-h_'+str(h)+'-b_'+str(b)+'-'+polarization+'-rotated'

"""### Output the dielectric profile with Ncellx x Ncelly unit cells 
### for zmin <= z <= zmax, with Nz values of z 
Ncellx = 5
Ncelly = 5
zmin = -0.2*Lz
zmax = 0.2*Lz 
Nx = 300
Ny = 300
Nz = 9

# Calculate the dielectric profile in planes parallel to Oxy
x_plot,y_plot,z_array,epsilon_xy_array \
    = DielectricProfileXY(ms,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz)
    
# Plot the dielectric profile, the name of the figure is:
#           namesave+'-z_'+str(k)+'.png'
# where k is the number of the value of z in the array z_array
PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave,'Yes')"""


### Plot the band structure
fig, ax = plt.subplots()
ax.plot(number, freqs)
plt.vlines(Nk+1,0.0,1.0,linestyle='dashed',color='black')
plt.vlines(2*(Nk+1),0.0,1.0,linestyle='dashed',color='black')
plt.xlim(0,3*(Nk+1))
plt.ylim(0,0.5)
tick_locs = [i*(Nk+1) for i in range(4)]
tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs,size=16)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
plt.title(namesave,fontsize=14)
plt.savefig(namesave+'.png')
plt.show()


### Print the band structure to file 
#PrintBandStructure(freqs,number,namesave)

with open(namesave+'-Band.txt','w') as file:
    for n in number:
        file.write('%.8f    ' % n)
        file.writelines('%.8f   ' % w for w in freqs[n])
        file.write('\n')