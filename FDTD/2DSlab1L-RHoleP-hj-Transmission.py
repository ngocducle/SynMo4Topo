import numpy as np
import scipy 
import meep as mp
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src')
from Materials import *
from FDTD_2Dstructures import geo_2DSlab1L_RholeP_hj_sameMater 

import os 

### Resolution 
resolution = 20 

### Boundary layers
# PML
dboundary = 3.0     # PML/Absorber thickness 
pml_layers = [mp.PML(direction = mp.X,
                     thickness = dboundary),
              mp.PML(direction = mp.Z,
                     thickness = dboundary)
            ]

abs_layers = [
    mp.Absorber(direction = mp.X,
                thickness = dboundary),
    mp.Absorber(direction = mp.Z,
                thickness = dboundary)
]

### The number of unit cells along the line y = 0 for each side
Ncell = 8

### Padding block
pad = 3.0 

### Geometrical parameters 
d = np.sqrt(2.0)    # The diagonal of one square unit cell 
h = 0.35            # Thickness of the layer 

# Left-hand side
b1 = 0.35            # The edge length of the undeformed square hole
e1 = -0.1             # The deformation parameter 

# Right-hand side
b2 = 0.41            # The edge length of the undeformed square hole
e2 = -0.1             # The deformation parameter 

# The total size of the structure along the x-axis
structurex = 2*Ncell*d 

# The toal size of the structure along the y-axis
structurey = d 

# The height of the unit cell along the z-direction (not including the PML/Absorber layer)
Lz = 3*h 

### Size of the simulation cell 
sx = structurex + 2*(pad + dboundary)
sy = d      # no PML, no pad, PBC 
sz = Lz + 2*dboundary

### Define the simulation cell 
cell = mp.Vector3(sx,sy,sz)

### The materials 
Mater = Si 
Envir = PMMA 

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

### The source 
fcen = 0.25    # pulse center frequency 
df = 0.05      # pulse width 
nfreq = 501    # number of frequencies
component = mp.Ey   # the component 
sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df), # is_integrated=False (default because no source extended into the PML)
        component = component,
        center = mp.Vector3(-0.5*sx+dboundary+0.5*pad,0,0),
        size = mp.Vector3(0,structurey,h)
    )
]

##### ===================================================================================
##### Compute the flux spectrum
freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0),
                     size = mp.Vector3(0,structurey,h),
                     direction = mp.X)

##### ====================================================================================
##### The position of the monitor
pt = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0)

##### =================================================================================
##### GEOMETRY 
geometry = geo_2DSlab1L_RholeP_hj_sameMater(d,h,
                                            vertice1,vertice2,
                                            Mater,Envir,Ncell,sx,sy)

##### ==================================================================================
##### Define the simulation 
sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = pml_layers,
        geometry = geometry,
        sources = sources,
        k_point = mp.Vector3(0,0,0), # PBC
	    ensure_periodicity = True,   # PBC
        resolution = resolution
    )

##### ==================================================================================
##### Transmitted flux
trans = sim.add_flux(fcen,df,nfreq,freg)

##### ==================================================================================
##### Run the simulation
sim.run(
    # Finish the simulation when the field at the point pt decays,
    # check the maximum of the field every dt
    until_after_sources = mp.stop_when_fields_decayed(
        dt = 50,
        c = component,
        pt = pt,
        decay_by = 1e-3
    )
)

##### ===================================================================================
### Get the dielectric function into array 
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
namesave = '2DSlab1L_'

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

##### ===============================================================================
##### Get the flux
trans_flux = np.array(mp.get_fluxes(trans))
freq_flux = np.array(mp.get_flux_freqs(trans))

datasave = np.column_stack((freq_flux,trans_flux))

print(np.shape(trans_flux))
print(np.shape(freq_flux))

##### ===============================================================================
##### Save the transmitted flux to file
filename = namesave + '.txt'
with open(filename,'w') as file:
    np.savetxt(file,datasave,'%.8f')