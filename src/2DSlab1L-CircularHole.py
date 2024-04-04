import numpy as np
import meep as mp 
from meep import mpb
import matplotlib.pyplot as plt

### The mode of the simulation, corresponding to the polarization
# Choose between: 'all', 'zeven', 'zodd' 
mode = 'zodd'

### Resolution 
resolution = mp.Vector3(16,16,16)   # pixels/a

### Geometrical parameters
h = 0.40     # Thickness of the slab
radius = 0.4 # Radius of the hole
Lz = 5.0     # The height of the unit cell along the z-direction

### Number of bands
num_bands = 18 

### Number of k-points to interpolate between 2 high-symmetry points
Nk = 9 

### Define the materials
Si = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0)
Air = mp.Medium(epsilon = 1.0)

Environment = Air 

### Define the lattice
geometry_lattice = mp.Lattice(
    size = mp.Vector3(1.0,1.0,Lz),
    basis1 = mp.Vector3(1.0,0.0),
    basis2 = mp.Vector3(0.0,1.0)
)

### Define the geometry
geometry = [
    mp.Block(
        center = mp.Vector3(0.0,0.0,0.0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Environment
    ),

    mp.Block(
        center = mp.Vector3(0.0,0.0,0.0),
        size = mp.Vector3(1.0,1.0,h),
        material = Si
    ),

    mp.Cylinder(
        center = mp.Vector3(0.0,0.0,0.0),
        radius = radius,
        height = h,
        axis = mp.Vector3(0,0,1),
        material = Environment
    )
]

### Define the k-points
# For the band structure in the whole Brillouin zone
k_points = [
    mp.Vector3(0.0,0.0,0.0),    # Gamma
    mp.Vector3(0.5,0.0,0.0),    # X
    mp.Vector3(0.5,0.5,0.0),    # M 
    mp.Vector3(0.0,0.0,0.0)     # Gamma 
]

# For the band in the vicinity of the M-point
#k_points = [
#    mp.Vector3(0.3,0.3,0.0),    # Gamma'
#    mp.Vector3(0.5,0.3,0.0),    # X
#    mp.Vector3(0.5,0.5,0.0),    # M 
#    mp.Vector3(0.3,0.3,0.0)     # Gamma 
#]

# Interpolate to get the points on the lines
k_points = mp.interpolate(Nk,k_points)

### Run the simulation
ms = mpb.ModeSolver(
    geometry = geometry,
    geometry_lattice = geometry_lattice,
    k_points = k_points,
    resolution = resolution,
    num_bands = num_bands
)

if mode == 'all':
    ms.run()
elif mode == 'zeven':
    ms.run_zeven()
elif mode == 'zodd':
    ms.run_zodd()
else:
    print('ERROR! The mode does not belong to the allowed list')
    exit()

freqs = ms.all_freqs

number = np.arange(len(k_points))

### Plot the bands 
fig, ax = plt.subplots()
ax.plot(number, freqs)
plt.show()