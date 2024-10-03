import numpy as np
import scipy 
import meep as mp 
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src')
from Heterojunctions import * 
from Materials import * 

import os 

### Resolution 
resolution = 20

### PML layer 
dpml = 2.0 # PML thickness 
pml_layers = [mp.PML(direction = mp.X,
                     thickness = dpml),
              mp.PML(direction = mp.Z,
                     thickness = dpml)    
            ]
#pml_layerZ = [mp.PML(direction = mp.Z,
#                     thickness = dpml)]

### The number of unit cells along the x-direction of each side 
Ncell = 3

### Padding block 
pad = 2.0 

### Geometrical parameters 
# The diagonal of one square unit cell
d = np.sqrt(2.0)

# The layer 1 
h1 = 0.80   # Thickness of the upper layer 
b1 = 0.32   # The edge length of the undeformed square hole 
e1 = -0.1   # The deformation parameter 

# The layer 2 
h2 = 0.40   # Thickness of the lower layer
b2 = 0.44   # The edge length of the undeformed square hole 
e2 = -0.1   # The deformation parameter 

# The interlayer distance 
dist = 1 

# The total size of the bilayer along the x-axis
structurex = (2*Ncell-0.5)*d 

# The total size of the bilayer along the y-axis
structurey = d 

# The total thickness of the bilayer 
hbilayer = h1 + dist + h2 

# The height of the unit cell along the z-direction
Lz = 1.5*hbilayer 

### Size of the simulation cell 
sx = (2*Ncell-0.5)*d + 2*(dpml+pad)
sy = d   # no PML, no pad 
sz = Lz + 2*dpml 

### Define the simulation cell 
cell = mp.Vector3(sx,sy,sz)

### The materials 
Mater = Si 
Envir = PMMA 
Mater1 = mp.Medium(epsilon=3.0) #Mater 
Mater2 = mp.Medium(epsilon=4.0) #Dielectric 

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

##### The array of shift 
### ATTENTION! Here is the fraction of the shift / d
Nq0 = 1
q0_array = np.linspace(0.0,0.0,Nq0)

### The source 
fcen = 0.258    # pulse center frequency  
df = 0.001      # pulse width 
nfreq = 501     # number of frequencies 
component = mp.Ey # the component 
sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df),
        component = component,
        center = mp.Vector3(-0.5*sx+dpml+0.5*pad,0,0),
        size = mp.Vector3(0,structurey,Lz)
    )
]

# The array of frequencies
freq_array = np.linspace(fcen-df,fcen+df,nfreq)

##### ==================================================================================
### The position of the monitor 
pt = mp.Vector3(0.5*sx-dpml-0.5*pad,0,0)

##### ==================================================================================
#####           We scan over the shift delta 
##### ==================================================================================
for idelta in range(Nq0):
    ### The shift 
    delta = d*q0_array[idelta]

    ##### ==================================================================================
    ##### GEOMETRY 
    geometry = [mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir
    )]

    geometry = geo_2DSlab2L_RHole_hj_PBCy(d,h1,b1,e1,h2,b2,e2,dist,
                                          vertice_cell,vertice1,vertice2,
                                          Mater1,Mater2,Envir,
                                          Ncell,sx,sy,
                                          structurex,structurey,hbilayer,delta)

    ##### ==================================================================================
    ##### Define the simulation 
    sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = pml_layers,
        geometry = geometry,
        sources = sources,
        resolution =resolution 
    )

    ##### ==================================================================================
    ##### Run the simulation 
    #sim.run(until_after_sources = mp.stop_when_fields_decayed(500,component,pt,1e-3))
    sim.run(until=50)

    ##### ==================================================================================
    ### Get the dielectric function into the array 
    eps_data = sim.get_array(center = mp.Vector3(0,0,0),
                             size = cell,
                             component = mp.Dielectric)

    print('sx = '+str(sx))
    print('sy = '+str(sy))
    print('sx = '+str(sz))
    print('Shape of eps_data: '+str(np.shape(eps_data)))

    ##### Plot the dielectric function 
    shape = np.shape(eps_data)
    Nx = shape[0]
    Ny = shape[1]
    Nz = shape[2]

    ##### The name of the files
    namesave = 'delta_{0:.4f}'.format(delta)

    os.system('mkdir '+namesave)

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

    os.system('mv *.png '+namesave)
