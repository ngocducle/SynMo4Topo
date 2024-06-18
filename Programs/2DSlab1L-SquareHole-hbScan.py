import numpy as np
import scipy  
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LSquareHole 
from BandStructure import * 
from ExportData import * 
from Materials import * 

##### The MAIN program goes here 
def main():
    #############################################################################
    #                                                                           #
    #       Please define the parameters of the simulation here                 #
    #                                                                           #
    #############################################################################


    ### The light polarization
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(8,8,8)  # pixels/a
    print('# The resolution:'+str(resolution))

    ### The array of thickness
    Nh = 5  
    h_array = np.linspace(0.05,0.50,Nh)

    ### The array of edge length of the square hole
    Nb = 5
    b_array = np.linspace(0.05,0.45,Nb)

    ### The height of the unit cell along the z-direction 
    Lz = 5.0 

    ### The number of numbands
    num_bands = 20 
    print('# The number of bands to simulate: '+str(num_bands))

    ### The set of k-points 
    Nk = 9          # The number of k-points to simulate 
    k_points = [
        mp.Vector3(0.4,0.4,0.0),
        mp.Vector3(0.5,0.4,0.0),
        mp.Vector3(0.5,0.5,0.0),
        mp.Vector3(0.4,0.4,0.0)
    ]

    k_points = mp.interpolate(Nk,k_points)

    ### Materials
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    #################################################################################
    #                                                                               #
    #                   Here the simulation starts                                  #
    #                                                                               #
    #################################################################################

    ##### We scan over the thickness and the length of the edge of the square hole
    for ih in range(Nh):
        for ib in range(Nb):
            h = h_array[ih]     # Thickness of the slab
            b = b_array[ib]     # Edge of the square hole 

            print('# *******************************************************')
            print('# Thickness of the slab h = '+str(h))
            print('# Edge of the square hole b = '+str(b))
            print('# The height of the simulation cell Lz = '+str(Lz))

            ### Define the mode solver 
            ms = _2DSlab1LSquareHole(h,Lz,b,num_bands,Nk,resolution,k_points,Mater,Envir)

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
            #namesave = '2DSlab1L-SquareHole-h_'+str(h)+'-b_'+str(b) \
            #    +'-'+polarization
            namesave = "2DSlab1L-SquareHole-h_{0:.2f}".format(h)
            namesave = namesave + "-b_{0:.2f}".format(b)
            namesave = namesave+'-'+polarization 
            print(namesave)

            ### Print the band structure to file 
            PrintBandStructure(freqs,number,namesave)

            ### Plot the band structure 
            PlotBand_M(number,freqs,Nk,namesave,'No')

##### Run the MAIN program 
if __name__ == "__main__":
    main()