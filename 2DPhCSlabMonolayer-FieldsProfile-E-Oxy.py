import numpy as np
import scipy 
import meep as mp
from meep import mpb 
import matplotlib.pyplot as plt 
from matplotlib import cm, colors 

### =======================================================================
### 
### Plot the field profiles in the Oxy plane 
### 
### =======================================================================

### Resolution
resolution = mp.Vector3(32,32,32) # pixels/a 

### Geometrical parameters 
h = 0.3      # Thickness of one layer 
radius = 0.4 # Radius of the hole
Lz = 5       # Length of the unit cell along the z-direction 

### The value of z where we take a slice to plot the field
zvalue = -0.25*h    
print('z = '+str(zvalue)) 

if abs(zvalue) > 0.5*Lz:
    print('ERROR! zvalue must be in the interval abs(zvalue) <= 0.5Lz')
    exit() 

### Number of bands 
num_bands = 20 

### Define the materials 
Si = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
PMMA = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0)
Air = mp.Medium(epsilon = 1.0)

Environment = PMMA 

### Define the lattice
geometry_lattice = mp.Lattice(size = mp.Vector3(1.0,1.0,Lz),  
                              basis1 = mp.Vector3(1.0,0.0),
                              basis2 = mp.Vector3(0.0,1.0))

### Define the geometry 
geometry = [ mp.Block(center = mp.Vector3(0.0,0.0,0.0), 
                      size = mp.Vector3(mp.inf,mp.inf,mp.inf), 
                      material = Environment),

             mp.Block(center = mp.Vector3(0.0,0.0,0.0), 
                      size = mp.Vector3(mp.inf,mp.inf,h), 
                      material = Si),  

             mp.Cylinder(center = mp.Vector3(0.0,0.0,0.0), 
                         radius = radius, 
                         height = h, 
                         axis = mp.Vector3(0.0,0.0,1.0), 
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

def get_fields(ms, band):
    efields.append(ms.get_efield(band,bloch_phase = True))

ms.run_zodd(
    mpb.output_at_kpoint(mp.Vector3(0.5,0.5,0.0),
                         mpb.fix_efield_phase,
                         get_fields
                        ) 
)

### Get the epsilon profile
resolution_eps = 81 # Number of pixels per a 
num_periods = 3 # Number of periods along each direction 
md = mpb.MPBData(rectify = True, 
                 periods = num_periods,  
                 resolution = resolution_eps)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)

print('The shape of converted_eps: '+str(np.shape(converted_eps)))

### ATTENTION! Be careful that the structure is also copied num_periods
# time along the z-axis 

### The epsilon profile in the Oxy plane (parallel to the plane Oxy: z = 0)
# Array of zvalues to plot eps
# The shape of converted_eps is 
# (num_periods*resolution_eps*length_unit_cell_x,
# num_periods*resolution_eps*length_unit_cell_y,
# num_periods*resolution_eps*length_unit_cell_z)
# Here:
# length_unit_cell_x = length_unit_cell_y = 1
# length_unit_cell_z = Lz 
shape_eps = np.shape(converted_eps)
len_zarray_eps = int(shape_eps[2] / num_periods)
print('len_zarray_eps = '+str(len_zarray_eps)) 
zarray_eps = np.linspace(-0.5*Lz,0.5*Lz,len_zarray_eps) 

# We look for the index of the element of zarray_eps which is closest
# to zvalue 
for i in range(len_zarray_eps-1): 
    if ((zarray_eps[i] <= zvalue) and (zvalue < zarray_eps[i+1])):
        zindex_eps = i 
        break 
    else:
        zindex_eps = len_zarray_eps - 1

# The index of the slice in the z-axis
print('zindex_eps = '+str(zindex_eps))

# The slice Oxy at z = 0
# The array eps_Oxy has shape
# (num_periods*resolution_eps*length_unit_cell_x,
# num_periods*resolution_eps*length_unit_cell_y) 
eps_Oxy = converted_eps[:,:,zindex_eps] 
print('The shape of eps_Oxy: '+str(np.shape(eps_Oxy))) 

# The numerical errors produce noise, the contours of the noise hide
# the field pattern, so we round eps_Oxy to 4 decimals 
eps_Oxy = np.round(eps_Oxy,4) 

# The meshgrid of (x,y)
Nx = resolution_eps*num_periods 
Ny = resolution_eps*num_periods 
Xlim = 0.5*num_periods   
Ylim = 0.5*num_periods  
X, Y = np.meshgrid( np.linspace(-Xlim,Xlim,Nx), 
                    np.linspace(-Ylim,Ylim,Ny) ) 

### Get the fields  
resolution_field = 81 
md = mpb.MPBData(rectify = True, 
                 resolution = resolution_field, 
                 periods = num_periods) 

# Define the arrays for X and Y to plot the fields
Nx = resolution_field*num_periods
Ny = resolution_field*num_periods

Xfield, Yfield = np.meshgrid( 
    np.linspace(-Xlim, Xlim, Nx),
    np.linspace(-Ylim, Ylim, Ny) 
    )

# Index for z: the z-coordinate of the slice 
# Array of the z-coordinate for field pattern (different from array of 
# dielectric constant)
zarray_field_len = int(Lz*resolution.z) 
zarray_field = np.linspace(-0.5*Lz,0.5*Lz,zarray_field_len)  

for i in range(zarray_field_len-1):
    if (zarray_field[i]<=zvalue) and (zvalue<zarray_field[i+1]):
        zindex_field = i 
        break
    else:
        zindex_field = zarray_field_len - 1 
  
print('zindex_field = '+str(zindex_field)) 

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
    # Get the x component of the E-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Ex at z
    Ex = f[:,:,zindex_field,0]  

    # Save the data Ex to Efieldx
    Efieldx.append(md.convert(Ex))  

    # Get the y component of the E-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Ey at z
    Ey = f[:,:,zindex_field,1]   

    # Save the data Ey to converted
    Efieldy.append(md.convert(Ey))  

    # Get the z component of the E-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Ez at z
    Ez = f[:,:,zindex_field,2] 

    # Save the data Ez to converted
    Efieldz.append(md.convert(Ez))  

for i in range(8):
    Ex = np.real(Efieldx[i])
    Ey = np.real(Efieldy[i])
    Ez = np.real(Efieldz[i]) 

    Xlim = int(num_periods/2)
    Ylim = int(num_periods/2) 

    fig, axs = plt.subplots(1,3)  
    axs[0].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[0].pcolormesh(Xfield,Yfield,Ex.T,shading='gouraud',cmap='RdBu')
    axs[0].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[0].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[0].set_xlabel('x/a',fontsize=14)
    axs[0].set_ylabel('y/a',fontsize=14) 
    axs[0].set_title('Ex', fontsize=14) 
    axs[0].set_aspect('equal') 

    axs[1].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[1].pcolormesh(Xfield,Yfield,Ey.T,shading='gouraud',cmap='RdBu') 
    axs[1].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[1].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[1].set_xlabel('x/a', fontsize=14)
    #axs[1].set_ylabel('y/a',fontsize=14) 
    axs[1].set_title('Ey', fontsize=14) 
    axs[1].set_aspect('equal') 

    axs[2].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[2].pcolormesh(Xfield,Yfield,Ez.T,shading='gouraud',cmap='RdBu')
    axs[2].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[2].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[2].set_xlabel('x/a', fontsize=14)
    #axs[2].set_ylabel('y/a', fontsize=14) 
    axs[2].set_title('Ez', fontsize=14) 
    axs[2].set_aspect('equal') 

    vmin = min(Ex.min(),Ey.min(),Ez.min()) 
    vmax = max(Ex.max(),Ey.max(),Ez.max())
    norm = colors.Normalize(vmin=vmin, vmax=vmax) 
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'),
                 orientation='vertical',  
                 shrink=0.4, 
                 ax=axs)      
 
    fig.suptitle('z = '+str(zvalue)+': Band '+str(i+1), fontsize=14)  
    plt.savefig('ReE_Oxy_Band'+str(i+1)+'.png')

#plt.show()   

for i in range(8):
    Ex = np.imag(Efieldx[i])
    Ey = np.imag(Efieldy[i])
    Ez = np.imag(Efieldz[i]) 

    Xlim = int(num_periods/2)
    Ylim = int(num_periods/2) 

    fig, axs = plt.subplots(1,3)  
    axs[0].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[0].pcolormesh(Xfield,Yfield,Ex.T,shading='gouraud',cmap='RdBu')
    axs[0].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[0].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[0].set_xlabel('x/a',fontsize=14)
    axs[0].set_ylabel('y/a',fontsize=14) 
    axs[0].set_title('Ex', fontsize=14) 
    axs[0].set_aspect('equal') 

    axs[1].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[1].pcolormesh(Xfield,Yfield,Ey.T,shading='gouraud',cmap='RdBu') 
    axs[1].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[1].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[1].set_xlabel('x/a', fontsize=14)
    #axs[1].set_ylabel('y/a',fontsize=14) 
    axs[1].set_title('Ey', fontsize=14) 
    axs[1].set_aspect('equal') 

    axs[2].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[2].pcolormesh(Xfield,Yfield,Ez.T,shading='gouraud',cmap='RdBu')
    axs[2].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[2].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[2].set_xlabel('x/a', fontsize=14)
    #axs[2].set_ylabel('y/a', fontsize=14) 
    axs[2].set_title('Ez', fontsize=14) 
    axs[2].set_aspect('equal') 

    vmin = min(Ex.min(),Ey.min(),Ez.min()) 
    vmax = max(Ex.max(),Ey.max(),Ez.max())
    norm = colors.Normalize(vmin=vmin, vmax=vmax) 
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'),
                 orientation='vertical',  
                 shrink=0.4, 
                 ax=axs)      
 
    fig.suptitle('z = '+str(zvalue)+': Band '+str(i+1), fontsize=14)  
    plt.savefig('ImE_Oxy_Band'+str(i+1)+'.png')

plt.show()   