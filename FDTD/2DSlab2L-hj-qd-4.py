import numpy as np
import scipy 
import meep as mp 
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src')
from Materials import * 

import os 

### Resolution 
resolution = 16

### PML layer 
dpml = 1.0  # PML thickness 
pml_layers = [mp.PML(dpml)]

### The number of unit cells along each direction 
# The number of unit cells along the horizontal direction for each half
Ncellx = 3

# The number of unit cells along the vertical direction 
Ncelly = 5 

# Padding block 
pad = 1.0 

### Geometrical parameters 
# The diagonal of one square unit cell 
d = np.sqrt(2.0) 

# The layer 1 
h1 = 0.35   # Thickness of the upper layer
b1 = 0.32   # The edge length of the undeformed square hole
e1 = -0.1   # The deformation parameter 

# The layer 2 
h2 = 0.35   # Thickness of the lower layer
b2 = 0.44   # The edge length of the undeformed square hole
e2 = -0.1   # The deformation parameter 

# The distant between the 2 layers 
dist = 0.1

# The total size of the bilayer along the x-axis
structurex = (2*Ncellx-0.5)*d 

# The total size of the bilayer along the y-axis
structurey = Ncelly*d 

# The total thickness of the bilayer 
htotal = h1 + dist + h2 

# The height of the unit cell along the z-direction 
Lz = 1.5*htotal

### The materials
Mater = Si 
Envir = PMMA 
Mater1 = Mater 
Mater2 = Mater 

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
    mp.Vector3(0,-b1*(1-e1)/(1+e1)/np.sqrt(2),0),
]

vertice2 = [
    mp.Vector3(b2*(1+e2)/(1-e2)/np.sqrt(2),0,0),
    mp.Vector3(0,b2*(1-e2)/(1+e2)/np.sqrt(2),0),
    mp.Vector3(-b2*(1+e2)/(1-e2)/np.sqrt(2),0,0),
    mp.Vector3(0,-b2*(1-e2)/(1+e2)/np.sqrt(2),0),
]

### Size of the simulation cell 
sx = (2*Ncellx - 0.5)*d + 2*(dpml+pad)
sy = Ncelly*d + 2*(dpml+pad)
sz = Lz + 2*dpml 

### Define the simulation cell 
cell = mp.Vector3(sx,sy,sz)

### The shift along the main diagonal 
delta = 0.0 

##### ==============================================================================
### The source 
fcen = 0.30         # pulse center frequency 
df   = 0.15         # pulse width 
nfreq = 101         # number of frequencies

sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df),
        component = mp.Ez,
        center = mp.Vector3(-0.5*sx+dpml+0.5*pad,0,0),
        size = mp.Vector3(0,structurey,Lz)
    )
]

# The array of frequencies 
freq_array = np.linspace(fcen-df,fcen+df,nfreq)

##### =============================================================================
##### Compute the flux spectrum 
freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dpml-0.5*pad),
                     size = mp.Vector3(0,2*structurey,Lz),
                     direction = mp.X)

### Simulation volum 
vol = mp.Volume(center = mp.Vector3(0,0,0),
                size = mp.Vector3(sx,sy,Lz))

### The position of the monitor 
pt = mp.Vector3(0.5*sx-dpml-0.5*pad,0,0)

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
        center = mp.Vector3(-0.25*sx,0,0.5*(htotal-h1)),
        size = mp.Vector3(0.5*sx,sy,h1),
        material = Mater1 
    )
)

geometry.append(
    mp.Block(
        center = mp.Vector3(-0.25*sx,0,0.5*(-htotal+h2)),
        size = mp.Vector3(0.5*sx,sy,h2),
        material = Mater2 
    )
)

geometry.append(
    mp.Block(
        center = mp.Vector3(0.25*sx,0,0.5*(htotal-h2)),
        size = mp.Vector3(0.5*sx,sy,h2),
        material = Mater2 
    )
)

geometry.append(
    mp.Block(
        center = mp.Vector3(0.25*sx,0,0.5*(-htotal+h1)),
        size = mp.Vector3(0.5*sx,sy,h1),
        material = Mater1 
    )
)

##### Add the unit cells at the heterojunction boundary 
for j in range(Ncelly):
    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.25*d,0.5*structurey-0.5*d-j*d,0.5*(htotal-h1)),
            material = Mater1 
        )
    )

    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.25*d,0.5*structurey-0.5*d-j*d,0.5*(-htotal+h2)),
            material = Mater2 
        )
    )

for j in range(Ncelly+1):
    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d,0.5*structurey-j*d,0.5*(htotal-h2)),
            material = Mater2
        )
    )

    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d,0.5*structurey-j*d,0.5*(-htotal+h1)),
            material = Mater1
        )
    )

##### Add the holes 
for j in range(Ncelly):
    for i in range(Ncellx):
        geometry.append(mp.Prism(
            vertices = vertice1,
            height=h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.5*structurex+0.5*d+i*d+0.5*delta,
                      0.5*structurey-0.5*d-j*d,
                      0.5*(htotal-h1)),
            material = Envir 
        ))

        geometry.append(mp.Prism(
            vertices = vertice2,
            height=h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.5*structurex+0.5*d+i*d-0.5*delta,
                      0.5*structurey-0.5*d-j*d,
                      0.5*(-htotal+h2)),
            material = Envir 
        ))

        
    for i in range(Ncellx):
        geometry.append(mp.Prism(
            vertices = vertice2,
            height=h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.75*d+i*d+0.5*delta,
                      0.5*structurey-0.5*d-j*d,
                      0.5*(htotal-h2)),
            material = Envir 
        ))

        geometry.append(mp.Prism(
            vertices = vertice1,
            height=h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.75*d+i*d-0.5*delta,
                      0.5*structurey-0.5*d-j*d,
                      0.5*(-htotal+h1)),
            material = Envir 
        ))

for j in range(Ncelly+1):
    for i in range(Ncellx):
        geometry.append(mp.Prism(
            vertices = vertice1,
            height=h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.5*structurex+i*d+0.5*delta,
                      0.5*structurey-j*d,
                      0.5*(htotal-h1)),
            material = Envir 
        ))

        geometry.append(mp.Prism(
            vertices = vertice2,
            height=h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.5*structurex+i*d-0.5*delta,
                      0.5*structurey-j*d,
                      0.5*(-htotal+h2)),
            material = Envir 
        ))

    for i in range(Ncellx+1):
        geometry.append(mp.Prism(
            vertices = vertice2,
            height=h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d+i*d+0.5*delta,
                      0.5*structurey-j*d,
                      0.5*(htotal-h2)),
            material = Envir 
        ))

        geometry.append(mp.Prism(
            vertices = vertice1,
            height=h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d+i*d-0.5*delta,
                      0.5*structurey-j*d,
                      0.5*(-htotal+h1)),
            material = Envir 
        ))

##### ===============================================================================
##### Define the simulation 
sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = pml_layers,
        geometry = geometry,
        sources = sources, 
        resolution = resolution
    )

##### Transmitted flux
trans = sim.add_flux(fcen,df,nfreq,freg)

##### ===============================================================================
##### Run the simulation 
sim.run(until_after_sources = mp.stop_when_fields_decayed(500,
                                                          mp.Ez,
                                                          pt,
                                                          1e-4))

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

os.system('mkdir structure')

for i in range(Nx):
    plt.figure()
    plt.imshow(eps_data[i,:,:].transpose(),interpolation='spline36',cmap='PuOr')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.savefig('x-'+str(i)+'.png')
    plt.close()

for j in range(Ny):
    plt.figure()
    plt.imshow(eps_data[:,j,:].transpose(),interpolation='spline36',cmap='PuOr')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.savefig('y-'+str(j)+'.png')
    plt.close()

for k in range(Nz):
    plt.figure()
    plt.imshow(eps_data[:,:,k].transpose(),interpolation='spline36',cmap='PuOr')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('z-'+str(k)+'.png')
    plt.close()

os.system('mv *.png structure')

##### ===============================================================================
##### Get the flux
trans_flux = np.array(mp.get_fluxes(trans))

datasave = np.column_stack((freq_array,trans_flux))

print(np.shape(trans_flux))

##### Save the transmitted flux to file
with open('transmission_flux.txt','w') as file:
    np.savetxt(file,datasave,'%.8f') 

