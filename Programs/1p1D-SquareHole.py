import numpy as np 
import scipy 
import meep as mp 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LSquareHole 
from DielectricProfile import * 
from BandStructure import * 
from Materials import * 

##### The MAIN program goes here
def main():
    ##########################################################################
    #                                                                        #
    #        Please define the parameters of the simulation here             #
    #                                                                        #
    ##########################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'all'
    print('# Polarization = '+polarization)

    ### The part of the momentum space to plot the band structure
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'Gamma-M' 
    print('# The part of the momentum space to simulate:'+kSpace) 

    ### Resolution 
    resolution = mp.Vector3(16,16,16)   # pixels/a 
    print('# The resolution:'+str(resolution)) 

    ### Geometrical parameters 
    h = 0.35    # Thickness of the slab
    edge = 0.4  # Edge of the square hole
    dist = 0.05 # Distance between the two layers 
    Lz = 10.0   # The height of the unit cell along the z-direction 

    print('# Thickness of the slab h = '+str(h))
    print('# Edge of the square hole b = '+str(edge))
    print('# The height of the simulation cell Lz = '+str(Lz))
    print('# Interlayer distance dist = '+str(dist))

    ### The array of shift delta 
    Ndelta = 51 
    delta_array = np.linspace(0.0,0.5,Ndelta)
    print('# The array of delta is: ')
    print(delta_array)

    ### Number of bands 
    num_bands = 25 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Number of k-points to interpolate between 2 high-symmetry points 
    Nk = 49 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### Show figure (Yes/No)
    show_fig = 'Yes'
    print('# Show the figure: '+str(show_fig))

    ### Materials 
    # The material forming the structure
    Mater = Si 

    # The environment
    Envir = PMMA 

    ##########################################################################
    #                                                                        #
    #                   Here the simulation starts                           #
    #                                                                        #
    ##########################################################################

    ### We run over the array of delta 
    for i in range(Ndelta):
        ### Extract the value of delta 
        delta = delta_array[i] 

        ### Define the mode solver 
        ms = _2DSlab2LSquareHole(h,Lz,edge,dist,delta,delta,
                                 num_bands,Nk,resolution,kSpace,Mater,Envir)
        
        ### Run the simulation 
        if polarization == 'all':
            ms.run()
        elif polarization == 'zeven':
            ms.run_zeven()
        elif polarization == 'zodd':
            ms.run_zodd()
        else:
            print('ERROR! The polarization does not belong to the allowed list')
            exit()

        ### Extract the frequencies of the modes from the ModeSolver
        freqs = ms.all_freqs 

        ### The number of elements in k_points 
        number = np.arange(len(ms.k_points))

        ### The title and the name of the files
        namesave = '1p1D-SquareHole-h_'+str(h)+'-b_'+str(edge) \
            +'-'+polarization+'-'+kSpace+'-delta_'+str(i)
        
        ### Print the band structure to file 
        with open(namesave+'-Band.txt','w') as file:
            for n in number:
                file.write('%.8f    ' % n)
                file.writelines('%.8f   ' % w for w in freqs[n])
                file.write('\n')
        
##### Run the MAIN program 
if __name__ == "__main__":
    main()