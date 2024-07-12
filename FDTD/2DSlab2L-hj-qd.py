import meep as mp 
import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src')
from Materials import * 

### Resolution
resolution = 16

### PML layer 
dpml = 1.0  # PML thickness 
pml_layers = [mp.PML(dpml)]

### The number unit cells along one direction 
N = 5 

### Padding block 
pad = 1.0 

### Geometrical parameters
# The layer 1
h1 = 0.35   # Thickness of the upper layer 
b1 = 0.38   # The edge length of the undeformed square hole 
e1 = 0.1    # The anisotropy between the two diagonals of the layer 1 

# The layer 2 
h2 = 0.35   # Thickness of the lower layer 
b2 = 0.38   # The edge length of the undeformed square hole 
e2 = 0.1    # The anisotropy between the two diagonals of the layer 2 

# The distant between the 2 layers 
dist = 0.1 

# The total thickness of the bilayer 
htotal = h1 + dist + h2 

# The height of the unit cell along the z-direction 
Lz = 3.0 

# The shifts along the x and y direction
deltax = 0.1
deltay = 0.1 

### The materials 
Mater = Si 
Envir = PMMA 

### The vertices of the rhombus holes 
vertices1 = [
    mp.Vector3( 0.5*b1*(1+e1)/(1-e1)+0.5*deltax,  
                0.5*b1*(1+e1)/(1-e1)+0.5*deltay, 0.0 ),
    mp.Vector3( 0.5*b1*(1-e1)/(1+e1)+0.5*deltax, 
               -0.5*b1*(1-e1)/(1+e1)+0.5*deltay, 0.0 ),
    mp.Vector3(-0.5*b1*(1+e1)/(1-e1)+0.5*deltax, 
               -0.5*b1*(1+e1)/(1-e1)+0.5*deltay, 0.0 ),
    mp.Vector3(-0.5*b1*(1-e1)/(1+e1)+0.5*deltax,  
                0.5*b1*(1-e1)/(1+e1)+0.5*deltay, 0.0 ) 
]

vertices2 = [
    mp.Vector3( 0.5*b2*(1+e2)/(1-e2)-0.5*deltax,
        0.5*b2*(1+e2)/(1-e2)-0.5*deltay, 0.0 ),
    mp.Vector3( 0.5*b2*(1-e2)/(1+e2)-0.5*deltax, 
        -0.5*b2*(1-e2)/(1+e2)-0.5*deltay, 0.0 ),
    mp.Vector3(-0.5*b2*(1+e2)/(1-e2)-0.5*deltax, 
        -0.5*b2*(1+e2)/(1-e2)-0.5*deltay, 0.0 ),
    mp.Vector3(-0.5*b2*(1-e2)/(1+e2)-0.5*deltax,  
        0.5*b2*(1-e2)/(1+e2)-0.5*deltay, 0.0 ) 
]

### Size of the simulation cell 
sx = N + 2*dpml + 2*pad 
sy = N + 2*dpml + 2*pad 
sz = h1 + h2 + dist + 2*pad + 2*dpml 

### Define the simulation cell
cell = mp.Vector3(sx,sy,sz)

### Define the geometry 
geometry = [] 

### The position of the center of the starting unit cell from left and below
if N % 2 == 1: 
    Nstart = -(N-1)/2
else: 
    Nstart = -N/2+0.5

### Initialize the geometry with environment
geometry.append(mp.Block(
    center = mp.Vector3(0.0,0.0,0.0),
    size = mp.Vector3(mp.inf,mp.inf,mp.inf),
    material = Envir 
))

Mater1 = mp.Medium(epsilon = 5.0)
Mater2 = Si 

for i in range(N):
    ### Lower left 
    for j in np.arange(0,i+1): 
        # Upper layer 
        geometry.append(mp.Block(
            center = mp.Vector3(Nstart+i,
                                Nstart+j,
                                0.5*(htotal-h1)),
            size = mp.Vector3(1,1,h1),
            material = Mater2 
        ))

        geometry.append(mp.Prism(
            vertices = vertices1,
            height = h1,
            axis = mp.Vector3(0.0,0.0,1.0),
            center = (Nstart+i+0.5*deltax,
                      Nstart+j+0.5*deltay,
                      0.5*(htotal-h1)),
            material = Envir 
        ))

        # Lower layer 
        geometry.append(mp.Block(
            center = mp.Vector3(Nstart+i,
                                Nstart+j,
                                -0.5*(htotal-h2)),
            size = mp.Vector3(1,1,h2),
            material = Mater1 
        ))

        geometry.append(mp.Prism(
            vertices = vertices2,
            height = h2,
            axis = mp.Vector3(0.0,0.0,1.0),
            center = (Nstart+i-0.5*deltax,Nstart+j-0.5*deltay,-0.5*(htotal-h2)),
            material = Envir
        ))

    ### Upper right
    for j in range(i+1,N):
        # Upper layer 
        geometry.append(mp.Block(
            center = mp.Vector3(Nstart+i,
                                Nstart+j,
                                0.5*(htotal-h2)),
            size = mp.Vector3(1,1,h2),
            material = Mater1 
        ))

        geometry.append(mp.Prism(
            vertices = vertices2,
            height = h2,
            axis = mp.Vector3(0.0,0.0,1.0),
            center = (Nstart+i+0.5*deltax,Nstart+j+0.5*deltay,0.5*(htotal-h2)),
            material = Envir
        ))

        # Lower layer 
        geometry.append(mp.Block(
            center = mp.Vector3(Nstart+i,
                                Nstart+j,
                                -0.5*(htotal-h1)),
            size = mp.Vector3(1,1,h1),
            material = Mater2 
        ))

        geometry.append(mp.Prism(
            vertices = vertices1,
            height = h1,
            axis = mp.Vector3(0.0,0.0,1.0),
            center = (Nstart+i-0.5*deltax,Nstart+j-0.5*deltay,-0.5*(htotal-h1)),
            material = Envir
        ))

### Padding blocks
# Left 
geometry.append(mp.Block(
    center = mp.Vector3(-0.5*sx+0.5*(pad+dpml),0,0.5*(htotal-h1)),
    size = mp.Vector3(pad+dpml,sy,h1),
    material = Mater1 
))

geometry.append(mp.Block(
    center = mp.Vector3(-0.5*sx+0.5*(pad+dpml),0,0.5*(-htotal+h2)),
    size = mp.Vector3(pad+dpml,sy,h2),
    material = Mater2  
))

# Right 
geometry.append(mp.Block(
    center = mp.Vector3(0.5*sx-0.5*(pad+dpml),0,0.5*(htotal-h2)),
    size = mp.Vector3(pad+dpml,sy,h2),
    material = Mater2  
))

geometry.append(mp.Block(
    center = mp.Vector3(0.5*sx-0.5*(pad+dpml),0,0.5*(-htotal+h1)),
    size = mp.Vector3(pad+dpml,sy,h1),
    material = Mater1 
))

# Below
geometry.append(mp.Block(
    center = mp.Vector3(0,0.5*sy-0.5*(pad+dpml),0.5*(htotal-h2)),
    size = mp.Vector3(N,pad+dpml,h2),
    material = Mater1 
))

geometry.append(mp.Block(
    center = mp.Vector3(0,0.5*sy-0.5*(pad+dpml),0.5*(-htotal+h1)),
    size = mp.Vector3(N,pad+dpml,h1),
    material = Mater2 
))

# Above 
geometry.append(mp.Block(
    center = mp.Vector3(0,-0.5*sy+0.5*(pad+dpml),0.5*(htotal-h1)),
    size = mp.Vector3(N,pad+dpml,h1),
    material = Mater2 
))

geometry.append(mp.Block(
    center = mp.Vector3(0,-0.5*sy+0.5*(pad+dpml),0.5*(-htotal+h2)),
    size = mp.Vector3(N,pad+dpml,h2),
    material = Mater1  
))

##### Define the simulation
sim = mp.Simulation(
    cell_size = cell,
    boundary_layers = pml_layers,
    geometry = geometry,
    resolution = resolution 
)

### Run the simulation
sim.run(until = 200)

### Get the dielectric function into an array
eps_data = sim.get_array(center = mp.Vector3(0,0,0),
                         size = cell,
                         component = mp.Dielectric)

print('sx = '+str(sx))
print('sy = '+str(sy))
print('sz = '+str(sz))
print('Shape of eps_data: '+str(np.shape(eps_data)))

### Plot the dielectric function 
shape = np.shape(eps_data)
Nx = shape[0]
Ny = shape[1]
Nz = shape[2]

for i in range(Nx):
    plt.figure()
    plt.imshow(eps_data[i,:,:].transpose(),interpolation='spline36',cmap='binary')
    plt.axis('off')
    plt.savefig('i-'+str(i)+'.png')
    plt.close()

for k in range(Nz):
    plt.figure()
    plt.imshow(eps_data[:,:,k].transpose(),interpolation='spline36',cmap='binary')
    plt.axis('off')
    plt.savefig('k-'+str(k)+'.png')
    plt.close()