import numpy as np
import meep as mp
from meep import mpb
import matplotlib.pyplot as plt

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to plot the dielectric profile          ###
###                                                                                ###
### ============================================================================== ###

def DielectricProfileZ(ModeSolver,x,y,zmin,zmax,Nz):
    z_array = np.linspace(zmin,zmax,Nz)

    epsilon_array = np.zeros(Nz)

    for i in range(Nz):
        epsilon_array[i] = ModeSolver.get_epsilon_point(mp.Vector3(x,y,z_array[i]))

    return z_array,epsilon_array

def PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave):
    fig,ax = plt.subplots()
    plt.plot(z_array,epsilon_z_array)
    plt.xlabel('z',fontsize=14)
    plt.ylabel(r'$\epsilon$', fontsize = 14)
    plt.title('x = '+str(x)+', y = '+str(y), fontsize = 14) 
    plt.savefig(namesave+'-epsilon-z.png')     
    plt.show() 


