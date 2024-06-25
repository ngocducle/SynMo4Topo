import numpy as np
import meep as mp
from meep import mpb
import matplotlib.pyplot as plt

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to plot the dielectric profile          ###
###                                                                                ###
### ============================================================================== ###

##### FUNCTION: Calculate the dielectric profile along the line (x,y,z)
###   where zmin <= z <= zmax, with Nz values of z
###   Give the array of z: z_array
###   and the corresponding array of dielectric constant: epsilon_z_array
def DielectricProfileZ(ModeSolver,x,y,zmin,zmax,Nz):
    z_array = np.linspace(zmin,zmax,Nz)

    epsilon_z_array = np.zeros(Nz)

    for i in range(Nz):
        epsilon_z_array[i] = ModeSolver.get_epsilon_point(mp.Vector3(x,y,z_array[i]))

    return z_array,epsilon_z_array

##### FUNCTION: Plot the dielectric profile along the line (x,y,z)
###   where zmin <= z <= zmax, with Nz values of z 
###   The figure name is: namesave+'-epsilon-z.png'
def PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave,show_fig):
    fig,ax = plt.subplots()
    plt.plot(z_array,epsilon_z_array)
    plt.xlabel('z',fontsize=14)
    plt.ylabel(r'$\epsilon$', fontsize = 14)
    plt.title('x = '+str(x)+', y = '+str(y), fontsize = 14) 
    plt.savefig(namesave+'-epsilon-z.png')     

    if show_fig == 'Yes':
        plt.show() 
    plt.close()

##### FUNCTION: Calculate the dielectric profile on the planes parallel to the Oxy plane
###   We extract the values of epsilon at the points (x,y,z) where:
###             -0.5*Ncellx <= x <= 0.5*Ncellx
###             -0.5*Ncelly <= y <= 0.5*Ncelly
###             zmin <= z <= zmax
###
###   The plots are centered at the point (0,0,z)
###   The size of the grids are (Nx,Ny,Nz)
def DielectricProfileXY(ModeSolver,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz):
    x_array = np.linspace(-0.5*Ncellx,0.5*Ncelly,Nx)
    y_array = np.linspace(-0.5*Ncelly,0.5*Ncelly,Ny)
    z_array = np.linspace(zmin,zmax,Nz)

    x_plot = np.zeros((Nx, Ny))
    for j in range(Ny):
        x_plot[:,j] = x_array  

    y_plot = np.zeros((Nx, Ny))
    for i in range(Nx):
        y_plot[i,:] = y_array 

    epsilon_xy_array = np.zeros((Nx,Ny,Nz))

    for k in range(Nz):
        for i in range(Nx):
            for j in range(Ny):
                epsilon_xy_array[i,j,k] = ModeSolver.get_epsilon_point \
                    (mp.Vector3(x_array[i],y_array[j],z_array[k]))
                     
    return x_plot,y_plot,z_array,epsilon_xy_array

##### FUNCTION: Plot the dielectric profile on the planes parallel to the Oxy plane
###   We extract the values of epsilon at the points (x,y,z) where:
###             -0.5*Ncellx <= x <= 0.5*Ncellx
###             -0.5*Ncelly <= y <= 0.5*Ncelly
###             zmin <= z <= zmax
###
###   The plots are centered at the point (0,0,z)
###   The size of the grids are (Nx,Ny,Nz)
def PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave,show_fig):
    Nz = len(z_array)
    vmin = 0 
    vmax = epsilon_xy_array.max()+0.1

    for k in range(Nz):
        fig,ax = plt.subplots()
        plt.pcolormesh(x_plot,y_plot,epsilon_xy_array[:,:,k],
                      cmap='RdBu',vmin=vmin,vmax=vmax)
        plt.xlabel('x',fontsize=14)
        plt.ylabel('y',fontsize=14)
        plt.title('z = '+str(z_array[k]),fontsize=14)
        plt.colorbar()
        ax.set_aspect('equal')
        plt.savefig(namesave+'-z_'+str(k)+'.png')
        
        if show_fig == 'Yes':
            plt.show()
        plt.close()

##### FUNCTION: Get the dielectric profile of the plane z = zvalue
###   Return the results as a MPBArray 
def DielectricProfileZvalue(ModeSolver,zvalue,Lz,resolution_eps,num_periods):
    ### Get the epsilon profile
    md = mpb.MPBData(rectify = True, 
                     periods = num_periods,  
                     resolution = resolution_eps)
    eps = ModeSolver.get_epsilon()
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

    # The slice Oxy at z = zvalue
    # The array eps_Oxy has shape
    # (num_periods*resolution_eps*length_unit_cell_x,
    # num_periods*resolution_eps*length_unit_cell_y) 
    eps_Oxy = converted_eps[:,:,zindex_eps] 
    print('The shape of eps_Oxy: '+str(np.shape(eps_Oxy)))

    # The numerical errors produce noise, the contours of the noise hide
    # the field pattern, so we round eps_Oxy to 4 decimals 
    eps_Oxy = np.round(eps_Oxy,4) 

    return eps_Oxy