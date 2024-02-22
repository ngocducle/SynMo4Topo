import numpy as np
import scipy 
import meep as mp
from meep import mpb 
import matplotlib.pyplot as plt

### =======================================================================
### 
### Plot the dielectric constant profile in the Oxz plane 
### 
### =======================================================================

### Resolution
resolution = mp.Vector3(32,32,32) # pixels/a

### Geometrical parameters 
h = 0.3      # Thickness of one layer 
radius = 0.4 # Radius of the hole
Lz = 5       # Length of the unit cell along the z-direction 

### Number of bands 
num_bands = 8 

### Define the materials 
Si = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
PMMA = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0)
Air = mp.Medium(epsilon = 1.0)

Environment = PMMA 

### Define the lattice
geometry_lattice = mp.Lattice(size = mp.Vector3(1, 1, Lz),  
                              basis1 = mp.Vector3(1, 0),
                              basis2 = mp.Vector3(0, 1))

### Define the geometry 
geometry = [ mp.Block(center = mp.Vector3(0, 0, 0), 
                      size = mp.Vector3(mp.inf, mp.inf, mp.inf), 
                      material = Environment),

             mp.Block(center = mp.Vector3(0, 0, 0), 
                      size = mp.Vector3(1, 1, h), 
                      material = Si),  

             mp.Cylinder(center = mp.Vector3(0, 0, 0), 
                         radius = radius, 
                         height = h, 
                         axis = mp.Vector3(0, 0, 1), 
                         material = Environment)
           ] 

# The k-point
k_points = [mp.Vector3(0.5,0.5,0.0)]

### Run the simulation 
ms = mpb.ModeSolver(
                    geometry = geometry,
                    geometry_lattice = geometry_lattice,
                    k_points = k_points,
                    resolution = resolution,
                    num_bands = num_bands 
)

### Run the z-even mode
ms.run_zeven(
    mpb.output_at_kpoint(mp.Vector3(0.5,0.5,0.0),
                         mpb.fix_hfield_phase 
                        ) 
)

### Get the epsilon profile
resolution_eps = 64 # Number of pixels per a 
num_periods = 1 # Number of periods along each direction 
md = mpb.MPBData(rectify = True, 
                 periods = num_periods,  
                 resolution = resolution_eps)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)

print('The shape of converted_eps: '+str(np.shape(converted_eps)))

### The epsilon profile in the Oyz plane (x = 0)
# The index of the slice in the x-axis
xindex = int(resolution_eps * num_periods / 2) 
print('zindex = '+str(xindex))

# The slice Oyz at x = 0
#eps_Oyz = np.swapaxes(converted_eps,1,2)
eps_Oyz = converted_eps[xindex,...] 
print('The shape of converted_eps: '+str(np.shape(converted_eps)))
print('The shape of eps_Oxz: '+str(np.shape(eps_Oyz)))

# The meshgrid of (y,z)
Ny = resolution_eps * num_periods 
Nz = Lz * resolution_eps * num_periods 
Ylim = 0.5 * num_periods   
Zlim = 0.5 * num_periods * Lz   
Y, Z = np.meshgrid( np.linspace(-Ylim,Ylim,Ny), np.linspace(-Zlim,Zlim,Nz) ) 

# Plot the epsilon-profile
# ATTENTION! Don't forget to transpose eps_Oxz 
fig, ax = plt.subplots()
plt.contour(Y,Z,eps_Oyz.T,cmap='binary')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14) 
ax.set_xlabel('y/a',fontsize=14)
ax.set_ylabel('z/a',fontsize=14) 
plt.title('x = 0', fontsize = 14) 
ax.set_aspect('equal') 
plt.show()

