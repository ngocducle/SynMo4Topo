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

##### FUNCTION: Print dielectric profile along z-direction at fixed (x,y)
def PrintDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave):
    Nz = len(z_array)

    with open(namesave+'-epsilon-z.txt','w') as file:
        for i in range(Nz):
            file.write('%.8f     ' % x)
            file.writelines('%.8f     ' % y)
            file.writelines('%.8f      ' % z_array[i])
            file.writelines('%.8f      ' % epsilon_z_array[i])
            file.write('\n') 