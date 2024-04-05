import numpy as np

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to print the data to files              ###
###                                                                                ###
### ============================================================================== ###

##### FUNCTION: Print band structure
def PrintBandStructure(freqs,number,namesave):
    with open(namesave+'-Band.txt','w') as file:
        for n in number:
            file.write('%.8f    ' % n)
            file.writelines('%.8f   ' % w for w in freqs[n])
            file.write('\n')

##### FUNCTION: Print the dielectric profile along the line (x,y,z)
###   where zmin <= z <= zmax, with Nz values of z 
###   The name of the data file is: namesave+'-epsilon-z.txt'
def PrintDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave):
    Nz = len(z_array)

    with open(namesave+'-epsilon-z.txt','w') as file:
        for k in range(Nz):
            file.write('%.8f     ' % x)
            file.writelines('%.8f     ' % y)
            file.writelines('%.8f      ' % z_array[k])
            file.writelines('%.8f      ' % epsilon_z_array[k])
            file.write('\n') 

##### FUNCTION: Print the dielectric profile on the planes parallel to the Oxy plane
###   We extract the values of epsilon at the points (x,y,z) where:
###             -0.5*Ncellx <= x <= 0.5*Ncellx
###             -0.5*Ncelly <= y <= 0.5*Ncelly
###             zmin <= z <= zmax
###
###   The plots are centered at the point (0,0,z)
###   The size of the grids are (Nx,Ny,Nz)
def PrintDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave):
    (Nx,Ny) = np.shape(x_plot)
    
    Nz = len(z_array)

    with open(namesave+'-epsilon-xy.txt','w') as file:
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    file.write('%.8f    ' % x_plot[i,j])
                    file.writelines('%.8f   ' % y_plot[i,j])
                    file.writelines('%.8f   ' % z_array[k])
                    file.writelines('%.8f   ' % epsilon_xy_array[i,j,k])
                    file.write('\n')

