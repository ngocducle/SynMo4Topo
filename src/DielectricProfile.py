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
def PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave):
    fig,ax = plt.subplots()
    plt.plot(z_array,epsilon_z_array)
    plt.xlabel('z',fontsize=14)
    plt.ylabel(r'$\epsilon$', fontsize = 14)
    plt.title('x = '+str(x)+', y = '+str(y), fontsize = 14) 
    plt.savefig(namesave+'-epsilon-z.png')     
    plt.show() 

##### FUNCTION: Plot the dielectric profile on the planes parallel to the Oxy plane
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

    xplot = np.zeros((Nx, Ny))
    for j in range(Ny):
        xplot[:,j] = x_array  

    yplot = np.zeros((Nx, Ny))
    for i in range(Nx):
        yplot[i,:] = y_array 

    epsilon_xy_array = np.zeros((Nx,Ny,Nz))

    for k in range(Nz):
        for i in range(Nx):
            for j in range(Ny):
                epsilon_xy_array[i,j,k] = ModeSolver.get_epsilon_point \
                    (mp.Vector3(x_array[i],y_array[j],z_array[k]))
                     
    return x_array,y_array,z_array,epsilon_xy_array





