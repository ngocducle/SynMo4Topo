import numpy as np
import scipy 
import meep as mp
from meep import mpb 
import matplotlib.pyplot as plt
from matplotlib import cm, colors 

### =======================================================================
### 
### Plot the field profiles in the Oyz plane 
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

### Run the simulation  
efields = []
hfields = []

def get_fields(ms, band):
    efields.append(ms.get_efield(band,bloch_phase = True))
    hfields.append(ms.get_hfield(band,bloch_phase = True)) 

ms.run_zeven(
    mpb.output_at_kpoint(mp.Vector3(0.5,0.5,0.0),
                         mpb.fix_efield_phase,
                         mpb.fix_hfield_phase,
                         get_fields
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
print('xindex = '+str(xindex))

# The slice Oyz at x = 0
#eps_Oyz = np.swapaxes(converted_eps,1,2)
eps_Oyz = converted_eps[xindex,:,:] 
print('The shape of eps_Oyz: '+str(np.shape(eps_Oyz)))

# The meshgrid of (y,z)
Ny = resolution_eps * num_periods 
Nz = Lz * resolution_eps * num_periods 
Ylim = 0.5 * num_periods   
Zlim = 0.5 * num_periods * Lz   
Y, Z = np.meshgrid( np.linspace(-Ylim,Ylim,Ny), np.linspace(-Zlim,Zlim,Nz) ) 

### Get the fields  
resolution_field = 64 
md = mpb.MPBData(rectify = True, 
                 resolution = resolution_field, 
                 periods = num_periods) 

# efields and hfields contain num_bands array (datasheets?), 
# each corresponds to one band. 
# Each array has 4 indices (x,y,z,E-mp.Vector3)
# First, swap the axes y and z, i.e. swapping indices 1 and 2 
# We need to swap to move the indices of the y-coordinate and 
# of the field to the two last positions, so that the other indices
# are at the leftmost side to be able to be replaced by the ellipsis.
# This is because we can put only 1 ellipsis into the formula. 

### ******* PLOT THE E-fields *************

Efieldx = []
Efieldy = []
Efieldz = [] 

for f in efields:
    # Index for x = 0
    xindex = int(resolution.x / 2)
    print('xindex = '+str(xindex))

    # Get the x component of the E-fields
    print('The shape of f: '+str(np.shape(f)))

    # Take the slice (band,y,z) for Ex at x = 0
    Ex = f[xindex,:,:,0]

    # Save the data Ex to converted
    Efieldx.append(md.convert(Ex)) 

    # Get the y component of the E-fields
    print('The shape of f: '+str(np.shape(f)))

    # Take the slice (band,y,z) for Ey at x = 0
    Ey = f[xindex,:,:,1]

    # Save the data Ey to converted
    Efieldy.append(md.convert(Ey))

    # Get the z component of the E-fields
    print('The shape of f: '+str(np.shape(f)))

    # Take the slice (band,y,z) for Ez at x = 0
    Ez = f[xindex,:,:,2]

    # Save the data Ex to converted
    Efieldz.append(md.convert(Ez))   

# Define the arrays for Y and Z to plot the fields
Yfield, Zfield = np.meshgrid( np.linspace(-0.5,0.5,resolution_field),
                              np.linspace(-0.5*Lz,0.5*Lz,resolution_field) )

for i in range(8):
    Ex = np.real(Efieldx[i])
    Ey = np.real(Efieldy[i])
    Ez = np.real(Efieldz[i]) 

    Zlim = int(Lz/2) 

    fig, axs = plt.subplots(1,3)  
    axs[0].contour(Y,Z,eps_Oyz.T,cmap='binary') 
    axs[0].pcolormesh(Yfield,Zfield,Ex.T,shading='gouraud',cmap='RdBu')
    axs[0].set_xticks(np.linspace(-0.5,0.5,2)) 
    axs[0].set_yticks(np.linspace(-Zlim,Zlim,int(Lz)+1)) 
    axs[0].set_xlabel('y/a',fontsize=14)
    axs[0].set_ylabel('z/a',fontsize=14) 
    axs[0].set_title('Ex', fontsize=14) 
    axs[0].set_aspect('equal') 

    axs[1].contour(Y,Z,eps_Oyz.T,cmap='binary') 
    axs[1].pcolormesh(Yfield,Zfield,Ey.T,shading='gouraud',cmap='RdBu') 
    axs[1].set_xticks(np.linspace(-0.5,0.5,2)) 
    axs[1].set_yticks(np.linspace(-Zlim,Zlim,int(Lz)+1)) 
    axs[1].set_xlabel('y/a', fontsize=14)
    #axs[1].set_ylabel('z/a',fontsize=14) 
    axs[1].set_title('Ey', fontsize=14) 
    axs[1].set_aspect('equal') 

    axs[2].contour(Y,Z,eps_Oyz.T,cmap='binary') 
    axs[2].pcolormesh(Yfield,Zfield,Ez.T,shading='gouraud',cmap='RdBu')
    axs[2].set_xticks(np.linspace(-0.5,0.5,2)) 
    axs[2].set_yticks(np.linspace(-Zlim,Zlim,int(Lz)+1)) 
    axs[2].set_xlabel('y/a', fontsize=14)
    #axs[2].set_ylabel('z/a', fontsize=14) 
    axs[2].set_title('Ez', fontsize=14) 
    axs[2].set_aspect('equal') 

    vmax = max(abs(Ex).max(),abs(Ey).max(),abs(Ez).max())
    norm = colors.Normalize(vmin=-vmax, vmax=vmax) 
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'),
                 orientation='vertical',  
                 shrink=1.0,  
                 ax=axs)      
 
    fig.suptitle('Band '+str(i+1)+': x = 0', fontsize=14)  
    plt.savefig('E_Oyz_Band'+str(i+1)+'.png')

plt.show()  



### ******* PLOT THE H-fields *************

Hfieldx = []
Hfieldy = []
Hfieldz = [] 

for f in hfields:
    # Index for y = 0
    xindex = int(resolution.y / 2)
    print('xindex = '+str(xindex))

    # Get the x component of the H-fields
    print('The shape of f: '+str(np.shape(f)))

    # Take the slice (band,y,z) for Hx at x = 0
    Hx = f[xindex,:,:,0]

    # Save the data Hx to converted
    Hfieldx.append(md.convert(Hx)) 

    # Get the y component of the H-fields
    print('The shape of f: '+str(np.shape(f)))

    # Take the slice (band,y,z) for Hy at x = 0
    Hy = f[xindex,:,:,1]

    # Save the data Hy to converted
    Hfieldy.append(md.convert(Hy))

    # Get the z component of the H-fields
    print('The shape of f: '+str(np.shape(f)))

    # Take the slice (band,y,z) for Hz at x = 0
    Hz = f[xindex,:,:,2]

    # Save the data Hx to converted
    Hfieldz.append(md.convert(Hz))   

# Define the arrays for X and Z to plot the fields
Yfield, Zfield = np.meshgrid( np.linspace(-0.5,0.5,resolution_field),
                              np.linspace(-0.5*Lz,0.5*Lz,resolution_field) )

for i in range(8):
    Hx = np.real(Hfieldx[i])
    Hy = np.real(Hfieldy[i])
    Hz = np.real(Hfieldz[i]) 

    Zlim = int(Lz/2) 

    fig, axs = plt.subplots(1,3)  
    axs[0].contour(Y,Z,eps_Oyz.T,cmap='binary') 
    axs[0].pcolormesh(Yfield,Zfield,Hx.T,shading='gouraud',cmap='RdBu')
    axs[0].set_xticks(np.linspace(-0.5,0.5,2)) 
    axs[0].set_yticks(np.linspace(-Zlim,Zlim,int(Lz)+1)) 
    axs[0].set_xlabel('y/a',fontsize=14)
    axs[0].set_ylabel('z/a',fontsize=14) 
    axs[0].set_title('Hx', fontsize=14) 
    axs[0].set_aspect('equal') 

    axs[1].contour(Y,Z,eps_Oyz.T,cmap='binary') 
    axs[1].pcolormesh(Yfield,Zfield,Hy.T,shading='gouraud',cmap='RdBu') 
    axs[1].set_xticks(np.linspace(-0.5,0.5,2)) 
    axs[1].set_yticks(np.linspace(-Zlim,Zlim,int(Lz)+1)) 
    axs[1].set_xlabel('y/a', fontsize=14)
    #axs[1].set_ylabel('z/a',fontsize=14) 
    axs[1].set_title('Hy', fontsize=14) 
    axs[1].set_aspect('equal') 

    axs[2].contour(Y,Z,eps_Oyz.T,cmap='binary') 
    axs[2].pcolormesh(Yfield,Zfield,Hz.T,shading='gouraud',cmap='RdBu')
    axs[2].set_xticks(np.linspace(-0.5,0.5,2)) 
    axs[2].set_yticks(np.linspace(-Zlim,Zlim,int(Lz)+1)) 
    axs[2].set_xlabel('y/a', fontsize=14)
    #axs[2].set_ylabel('z/a', fontsize=14) 
    axs[2].set_title('Hz', fontsize=14) 
    axs[2].set_aspect('equal') 

    vmax = max(abs(Hx).max(),abs(Hy).max(),abs(Hz).max())
    norm = colors.Normalize(vmin=-vmax, vmax=vmax) 
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'),
                 orientation='vertical',  
                 shrink=1.0,  
                 ax=axs)      
 
    fig.suptitle('Band '+str(i+1)+': y = 0', fontsize=14)  
    plt.savefig('H_Oyz_Band'+str(i+1)+'.png')

plt.show()  