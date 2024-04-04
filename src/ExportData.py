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

