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

# The component of the excitation 
component = mp.Ey 

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

##### The array of synthetic momentum 
Nq = 1
q_array = np.linspace(0.305,0.305,Nq)

### The shift along the main diagonal 
delta = 0.0 

##### ==============================================================================
### The source 
fcen = 0.258         # pulse center frequency 
df   = 0.001         # pulse width 
nfreq = 501          # number of frequencies

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
#####       We scan over the synthetic momentum q 
##### ==============================================================================
for iq in range(Nq):
    ### The shift 
    delta = q_array[iq]
    
    ##### ==============================================================================
    ##### GEOMETRY 
    ### Initialize the geometry with environment 
    geometry = geo_2DSlab2L_RHole_hj(d,h1,b1,e1,h2,b2,e2,dist,
                          vertice_cell,vertice1,vertice2,
                          Mater1,Mater2,Envir,
                          Ncellx,Ncelly,sx,sy,
                          structurex,structurey,htotal,delta)

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
                                                          component,
                                                          pt,
                                                          1e-3))

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

    ##### The name of the files
    namesave = 'q_{0:.4f}'.format(delta)

    os.system('mkdir '+namesave)

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

    os.system('mv *.png '+namesave)

    ##### ===============================================================================
    ##### Get the flux
    trans_flux = np.array(mp.get_fluxes(trans))

    datasave = np.column_stack((freq_array,trans_flux))

    print(np.shape(trans_flux))

    ##### Save the transmitted flux to file
    filename = namesave+'.txt'
    with open(filename,'w') as file:
        np.savetxt(file,datasave,'%.8f') 