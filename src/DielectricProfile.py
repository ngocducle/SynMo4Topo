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


