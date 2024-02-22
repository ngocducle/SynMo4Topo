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
resolution_eps = 81 # Number of pixels per a 
num_periods = 3 # Number of periods along each direction 
md = mpb.MPBData(rectify = True, 
                 periods = num_periods,  
                 resolution = resolution_eps)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)

print('The shape of converted_eps: '+str(np.shape(converted_eps)))

### The epsilon profile in the Oxy plane (z = 0)
# The index of the slice in the z-axis
zindex = int(Lz * resolution_eps * num_periods / 2) 
print('zindex = '+str(zindex))

# The slice Oxy at z = 0
eps_Oxy = converted_eps[...,zindex]
print('The shape of eps_Oxy: '+str(np.shape(eps_Oxy)))

# The meshgrid of (x,y)
Nx = resolution_eps * num_periods 
Ny = resolution_eps * num_periods 
Xlim = 0.5 * num_periods   
Ylim = 0.5 * num_periods  
X, Y = np.meshgrid( np.linspace(-Xlim,Xlim,Nx), np.linspace(-Ylim,Ylim,Ny) ) 

### Get the fields  
resolution_field = 81
md = mpb.MPBData(rectify = True, 
                 resolution = resolution_field, 
                 periods = num_periods) 

# Define the arrays for X and Y to plot the fields
Xfield, Yfield = np.meshgrid( 
    np.linspace(-Xlim, Xlim, Nx),
    np.linspace(-Ylim, Ylim, Ny) 
    )

# Index for z: the z-coordinate of the slice  
#zindex = int(Lz * resolution.z / 2) # z in the vicinity of 0 
zindex = 45   
print('zindex = '+str(zindex)) 

#converted = []

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
    Ex = f[:,:,zindex,0]  

    # Save the data Ex to converted
    #converted.append(md.convert(Ex)) 
    Efieldx.append(md.convert(Ex))  

    # Get the y component of the E-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Ey at z
    Ey = f[:,:,zindex,1]   

    # Save the data Ey to converted
    #converted.append(md.convert(Ey))
    Efieldy.append(md.convert(Ey))  

    # Get the z component of the E-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Ez at z
    Ez = f[:,:,zindex,2] 

    # Save the data Hz to converted
    #converted.append(md.convert(Ez))
    Efieldz.append(md.convert(Ez))  


# Plot the field-profile
#for i, Ex in enumerate(converted):
#    fig, ax = plt.subplots() 
#    plt.contour(X,Y,eps_Oxy.T,cmap='binary') 
    #plt.imshow(np.real(Ex).T,interpolation='spline36',cmap='RdBu',alpha=0.9)
#    plt.pcolormesh(X,Y,np.real(Ex).T,shading='gouraud',cmap='RdBu')
#    plt.colorbar() 
#    plt.xticks(fontsize = 14)
#    plt.yticks(fontsize = 14) 
#    ax.set_xlabel('x/a',fontsize=14)
#    ax.set_ylabel('y/a',fontsize=14) 
#    plt.title('Ex Band '+str(i+1)+': z = 0', fontsize = 14) 
#    ax.set_aspect('equal') 

#for i, Ey in enumerate(converted):
#    fig, ax = plt.subplots() 
#    plt.contour(X,Y,eps_Oxy.T,cmap='binary') 
    #plt.imshow(np.real(Ey).T,interpolation='spline36',cmap='RdBu',alpha=0.9)
#    plt.pcolormesh(X,Y,np.real(Ey).T,shading='gouraud',cmap='RdBu')
#    plt.colorbar() 
#    plt.xticks(fontsize = 14)
#    plt.yticks(fontsize = 14) 
#    ax.set_xlabel('x/a',fontsize=14)
#    ax.set_ylabel('y/a',fontsize=14) 
#    plt.title('Ey Band '+str(i+1)+': z = 0', fontsize = 14) 
#    ax.set_aspect('equal') 

#for i, Ez in enumerate(converted): 
#    fig, ax = plt.subplots()   
#    plt.contour(X,Y,eps_Oxy.T,cmap='binary') 
    #plt.imshow(np.real(Ez).T,interpolation='spline36',cmap='RdBu',alpha=0.9)
#    plt.pcolormesh(X,Y,np.real(Ez).T,shading='gouraud',cmap='RdBu')
#    plt.colorbar() 
#    plt.xticks(fontsize = 14)
#    plt.yticks(fontsize = 14) 
#    ax.set_xlabel('x/a',fontsize=14)
#    ax.set_ylabel('y/a',fontsize=14) 
#    plt.title('Ez Band '+str(i+1)+': z = 0', fontsize = 14) 
#    ax.set_aspect('equal') 

#    plt.show()

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
 
    fig.suptitle('Band '+str(i+1), fontsize=14)  
    plt.savefig('E_Oxy_Band'+str(i+1)+'.png')

#plt.show()   


### ******* PLOT THE H-fields *************

Hfieldx = []
Hfieldy = []
Hfieldz = [] 

for f in hfields:
    # Get the x component of the H-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Hx at z
    Hx = f[:,:,zindex,0]  

    # Save the data Hx to converted
    #converted.append(md.convert(Hx)) 
    Hfieldx.append(md.convert(Hx))  

    # Get the y component of the H-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Hy at z
    Hy = f[:,:,zindex,1]   

    # Save the data Hy to converted
    #converted.append(md.convert(Hy))
    Hfieldy.append(md.convert(Hy))  

    # Get the z component of the H-fields 
    print('The shape of f: '+str(np.shape(f))) 

    # Take the slice (band,x,y) for Hz at z
    Hz = f[:,:,zindex,2] 

    # Save the data Hz to converted
    #converted.append(md.convert(Hz))
    Hfieldz.append(md.convert(Hz))  


# Plot the field-profile
#for i, Hx in enumerate(converted):
#    fig, ax = plt.subplots() 
#    plt.contour(X,Y,eps_Oxy.T,cmap='binary') 
    #plt.imshow(np.real(Hx).T,interpolation='spline36',cmap='RdBu',alpha=0.9)
#    plt.pcolormesh(X,Y,np.real(Hx).T,shading='gouraud',cmap='RdBu')
#    plt.colorbar() 
#    plt.xticks(fontsize = 14)
#    plt.yticks(fontsize = 14) 
#    ax.set_xlabel('x/a',fontsize=14)
#    ax.set_ylabel('y/a',fontsize=14) 
#    plt.title('Hx Band '+str(i+1)+': z = 0', fontsize = 14) 
#    ax.set_aspect('equal') 

#for i, Hy in enumerate(converted):
#    fig, ax = plt.subplots() 
#    plt.contour(X,Y,eps_Oxy.T,cmap='binary') 
    #plt.imshow(np.real(Hy).T,interpolation='spline36',cmap='RdBu',alpha=0.9)
#    plt.pcolormesh(X,Y,np.real(Hy).T,shading='gouraud',cmap='RdBu')
#    plt.colorbar() 
#    plt.xticks(fontsize = 14)
#    plt.yticks(fontsize = 14) 
#    ax.set_xlabel('x/a',fontsize=14)
#    ax.set_ylabel('y/a',fontsize=14) 
#    plt.title('Hy Band '+str(i+1)+': z = 0', fontsize = 14) 
#    ax.set_aspect('equal') 

#for i, Hz in enumerate(converted): 
#    fig, ax = plt.subplots()   
#    plt.contour(X,Y,eps_Oxy.T,cmap='binary') 
    #plt.imshow(np.real(Hz).T,interpolation='spline36',cmap='RdBu',alpha=0.9)
#    plt.pcolormesh(X,Y,np.real(Hz).T,shading='gouraud',cmap='RdBu')
#    plt.colorbar() 
#    plt.xticks(fontsize = 14)
#    plt.yticks(fontsize = 14) 
#    ax.set_xlabel('x/a',fontsize=14)
#    ax.set_ylabel('y/a',fontsize=14) 
#    plt.title('Hz Band '+str(i+1)+': z = 0', fontsize = 14) 
#    ax.set_aspect('equal') 

#    plt.show()

for i in range(8):
    Hx = np.real(Hfieldx[i])
    Hy = np.real(Hfieldy[i])
    Hz = np.real(Hfieldz[i]) 

    Xlim = int(num_periods/2)
    Ylim = int(num_periods/2) 

    fig, axs = plt.subplots(1,3)  
    axs[0].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[0].pcolormesh(Xfield,Yfield,Hx.T,shading='gouraud',cmap='RdBu')
    axs[0].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[0].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[0].set_xlabel('x/a',fontsize=14)
    axs[0].set_ylabel('y/a',fontsize=14) 
    axs[0].set_title('Hx', fontsize=14) 
    axs[0].set_aspect('equal') 

    axs[1].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[1].pcolormesh(Xfield,Yfield,Hy.T,shading='gouraud',cmap='RdBu') 
    axs[1].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[1].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[1].set_xlabel('x/a', fontsize=14)
    #axs[1].set_ylabel('y/a',fontsize=14) 
    axs[1].set_title('Hy', fontsize=14) 
    axs[1].set_aspect('equal') 

    axs[2].contour(X,Y,eps_Oxy.T,cmap='binary') 
    axs[2].pcolormesh(Xfield,Yfield,Hz.T,shading='gouraud',cmap='RdBu')
    axs[2].set_xticks(np.linspace(-Xlim,Xlim,num_periods)) 
    axs[2].set_yticks(np.linspace(-Ylim,Ylim,num_periods)) 
    axs[2].set_xlabel('x/a', fontsize=14)
    #axs[2].set_ylabel('y/a', fontsize=14) 
    axs[2].set_title('Hz', fontsize=14) 
    axs[2].set_aspect('equal') 

    vmin = min(Hx.min(),Hy.min(),Hz.min()) 
    vmax = max(Hx.max(),Hy.max(),Hz.max())
    norm = colors.Normalize(vmin=vmin, vmax=vmax) 
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'),
                 orientation='vertical',  
                 shrink=0.4, 
                 ax=axs)      
 
    fig.suptitle('Band '+str(i+1), fontsize=14)  
    plt.savefig('H_Oxy_Band'+str(i+1)+'.png') 

plt.show()   