import numpy as np 
import scipy 
import meep as mp 
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src')
from Materials import * 

import os 

### Resolution 
resolution = 20

### PML layer 
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

### The number of unit cells along the line y = 0
Ncell = 8

### Padding block
pad = 3.0 

### Geometrical parameters 
# The diagonal of one square unit cell 
d = np.sqrt(2.0) 

# The layer 1 
h1 = 0.35   # Thickness of the upper layer
b1 = 0.41   # The edge length of the undeformed square hole
e1 = -0.1   # The deformation parameter 

# The layer 2 
h2 = 0.35   # Thickness of the lower layer
b2 = 0.35   # The edge length of the undeformed square hole
e2 = -0.1   # The deformation parameter 

# The distant between the 2 layers 
dist = 0.1

# The total thickness of the bilayer 
hbilayer = h1 + dist + h2 

# The height of the unit cell along the z-direction 
Lz = 3*hbilayer 

### Size of the simulation cell 
sx = 2*(pad+dboundary)
sy = d
sz = Lz + 2*dboundary 

### Define the simulation cell 
cell = mp.Vector3(sx,sy,sz)

### The materials
Mater = Si 
Envir = PMMA 

##### ====================================================================
### The source 
fcen = 0.26         # pulse center frequency 
df   = 0.01         # pulse width 
nfreq = 501          # number of frequencies
component = mp.Ey   # the component 

sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df), # is_integrated=False (default because no source extended into the PML)
        component = component,
        center = mp.Vector3(-0.5*sx+dboundary+0.5*pad,0,0),
        size = mp.Vector3(0,d,hbilayer)
    )
]

##### =============================================================================
##### Compute the flux spectrum 
freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dboundary-0.5*pad),
                     size = mp.Vector3(0,d,hbilayer),
                     direction = mp.X)

### Simulation volume
vol = mp.Volume(center = mp.Vector3(0,0,0),
                size = mp.Vector3(sx,sy,Lz))

### The position of the monitor 
pt = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0)

##### ==============================================================================
##### GEOMETRY
### Initialize the geometry with environment
geometry = []

geometry.append(mp.Block(
    center = mp.Vector3(0,0,0),
    size = mp.Vector3(mp.inf,mp.inf,mp.inf),
    material = Envir
))

### Add the structure (already includes the padding blocks)
geometry.append(
    mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,hbilayer),
        material = Mater 
    )
)

##### ===============================================================================
##### Define the simulation 
sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = pml_layers,
        geometry = geometry,
        sources = sources, 
        k_point = mp.Vector3(0,0,0), # PBC
        ensure_periodicity = True, # PBC
        resolution = resolution
    )

##### Transmitted flux
trans = sim.add_flux(fcen,df,nfreq,freg)

##### ===============================================================================
##### Run the simulation 
sim.run(
    # Finish the simulation when the field at the point pt decays,
    # check the maximum of the field every dt 
    until_after_sources = mp.stop_when_fields_decayed(
        dt = 3000,
        c = component,
        pt = pt,
        decay_by = 1e-4
    )
)

##### ====================================================================
### Get the dielectric function into an array 
eps_data = sim.get_array(center = mp.Vector3(0,0,0),
                         size = cell,
                         component = mp.Dielectric)

print('sx = '+str(sx))
print('sy = '+str(sy))
print('sz = '+str(sz))
print('Shape of eps_data: '+str(np.shape(eps_data)))

##### Plot the dielectric function
shape = np.shape(eps_data)
Nx = shape[0]
Ny = shape[1]
Nz = shape[2]

os.system('mkdir no_structure')

for i in range(Nx):
    plt.figure()
    plt.imshow(eps_data[i,:,:].transpose(),interpolation='spline36',cmap='coolwarm')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.savefig('x-'+str(i)+'.png')
    plt.close()

for j in range(Ny):
    plt.figure()
    plt.imshow(eps_data[:,j,:].transpose(),interpolation='spline36',cmap='coolwarm')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.savefig('y-'+str(j)+'.png')
    plt.close()

for k in range(Nz):
    plt.figure()
    plt.imshow(eps_data[:,:,k].transpose(),interpolation='spline36',cmap='coolwarm')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('z-'+str(k)+'.png')
    plt.close()

os.system('mv *.png no_structure')

##### ===============================================================================
##### Get the flux
trans_flux = np.array(mp.get_fluxes(trans))
freq_flux = np.array(mp.get_flux_freqs(trans))

datasave = np.column_stack((freq_flux,trans_flux))

print(np.shape(trans_flux))
print(np.shape(freq_flux))

##### Save the transmitted flux to file
with open('./no_structure/normalized_flux.txt','w') as file:
    np.savetxt(file,datasave,'%.8f') 