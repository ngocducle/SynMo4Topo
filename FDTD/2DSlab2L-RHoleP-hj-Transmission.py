import numpy as np
import scipy 
import meep as mp
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'../src')
from Materials import * 
from FDTD_2Dstructures import geo_2DSlab2L_RHoleP_hj_sameMater

import os 

### Resolution 
resolution = 20

### Boundary layers 
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

mix_layers = [mp.Absorber(direction = mp.X,
                          thickness = dboundary),
              mp.PML(direction = mp.Z,
                    thickness = dboundary)
            ]

### The source 
fcen = 0.260    # pulse center frequency 
df = 0.010      # pulse width 
nfreq = 1001     # number of frequencies

# The array of frequencies and wavelength
freq_array = np.linspace(fcen-df,fcen+df,nfreq)
wvl_array = 1/freq_array 

### The number of unit cells along the line y = 0
Ncell = 8

### Padding block
pad = 2.0*wvl_array[0] 

### Geometrical parameters 
d = np.sqrt(2.0)   # The diagonal of one square unit cell
h = 0.35           # Thickness of the two layers

# The layer 1 
b1 = 0.46   # The edge length of the undeformed square hole
e1 = -0.1    # The deformation parameter 

# The layer 2 
b2 = 0.30   # The edge length of the undeformed square hole 
e2 = -0.1    # The deformation parameter 

# The interlayer distance 
dist = 0.1

# The total size of the structure along the x-axis
structurex = (2*Ncell+1)*d 

# The toal size of the structure along the y-axis
structurey = d 

# The total thickness of the bilayer (z-direction)
hbilayer = 2*h + dist 

# The height of the unit cell along the z-direction (not including the PML/Absorber layer)
Lz = 3*hbilayer  

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

##### The array of shift 
### ATTENTION! Here is the fraction of the shift/d
Ndelta = 1 
delta_array = np.linspace(0.5,0.5,Ndelta)

##### The source
component = mp.Ey   # the component 
sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df,is_integrated=True),
        component = component,
        center = mp.Vector3(-0.5*sx+dboundary+0.5*pad,0,0),
        size = mp.Vector3(0,structurey,hbilayer)
    )
]

##### ===================================================================================
##### Compute the flux spectrum
freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0),
                     size = mp.Vector3(0,structurey,hbilayer),
                     direction = mp.X)

##### ====================================================================================
##### The position of the monitor
pt = mp.Vector3(0.5*sx-dboundary-0.5*pad,0,0)

##### =====================================================================================
#####               We scan over the shift 
##### =====================================================================================
for idelta in range(Ndelta):
    ### The shift 
    delta = d*delta_array[idelta]

    ##### ==================================================================================
    ##### GEOMETRY
    geometry = geo_2DSlab2L_RHoleP_hj_sameMater(d,h,hbilayer,delta,
                                                vertice1,vertice2,
                                                Mater,Envir,Ncell,sx,sy)

    ##### ==================================================================================
    ##### Define the simulation 
    sim = mp.Simulation(
        cell_size = cell,
        boundary_layers = abs_layers,
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
            dt = 100,
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
    namesave = 'delta_{0:.4f}'.format(delta_array[idelta])

    if mp.am_master():
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

    datasave = np.column_stack((freq_array,trans_flux))

    print(np.shape(trans_flux))

    ##### ===============================================================================
    ##### Save the transmitted flux to file
    if mp.am_master():
        filename = namesave + '.txt'
        with open(filename,'w') as file:
            np.savetxt(file,datasave,'%.8f')

