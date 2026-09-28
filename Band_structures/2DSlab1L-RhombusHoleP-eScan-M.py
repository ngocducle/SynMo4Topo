import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LRhombusHole
from Materials import * 

######################################################################################
#                                                                                    #
#       Calculation of the band structure for 2DSlab1L with rhombus hole             #
#       The hole is deformed from the C4 square, breaking the C4 symmetry            #
#       The product of the diagonal remains unchanged, so we add the suffix          #
#       P (product). It allows to keep the hole area, filling factor, and            #
#       effective refractive index unchanged                                         # 
#       That means let d be the diagonal of the square hole, the diagonals           #
#       of the rhombi are:                                                           #
#           d1 =  d*(1+e)/(1-e)                                                      #
#           d2 =  d*(1-e)/(1+e)                                                      #
#                                                                                    #
#       We fix the slab thickness h and the edge length b of the undeformed square   #
#       and scan the deformation parameter e                                         #
#       This program gives only the energy eigenvalue at the point M(0.5,0.5)        #
#                                                                                    # 
######################################################################################

##### The main program goes here 
def main():
    ##################################################################################
    #                                                                                #
    #          Please define the parameters of the simulation here                   #
    #                                                                                #
    ##################################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(16,16,16)   # pixels/a 

    ### Number of bands 
    num_bands = 20 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Geometrical parameters 
    h = 0.35        # Thickness of the slab
    Lz = 5.0        # The height of the unit cell along the z-direction 
    b = 0.4         # The mean of the projections of the diagonals vertices (0<=b<0.5)

    ### The array of anisotropy 
    Ne = 21
    e_array = np.linspace(-0.05,0.05,Ne) 

    print('# Thickness of the slab h = '+str(h))
    print('# The projection of mean half-axis of the rhombus on the edges b = '+str(b))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### Materials 
    # The material forming the structure
    Mater = Si 

    # The environment 
    Envir = PMMA 

    ### The set of k-points 
    k_points = [
        mp.Vector3(0.5,0.5,0.0) # M 
    ]

    # Show the figure 
    show_fig = 'Yes'

    #################################################################################
    #                                                                               #
    #                         Here the simulation starts                            #
    #                                                                               #
    #################################################################################

    ##### Scan in e_array 
    for i in range(Ne):
        ### The anisotropy 
        e = e_array[i]

        ### The vertices of the rhombus hole 
        # If e > 0 then the diagonal x=y is longer than the diagonal x=-y 
        # If e = 0 then the hole is a square 
        # If e < 0 then the diagonal x=y is shorter than the diagonal x=-y 
        vertices = [
            mp.Vector3(  0.5*b*(1+e)/(1-e),  0.5*b*(1+e)/(1-e), 0.0 ),
            mp.Vector3(  0.5*b*(1-e)/(1+e), -0.5*b*(1-e)/(1+e), 0.0 ),
            mp.Vector3( -0.5*b*(1+e)/(1-e), -0.5*b*(1+e)/(1-e), 0.0 ),
            mp.Vector3( -0.5*b*(1-e)/(1+e),  0.5*b*(1-e)/(1+e), 0.0 )
        ]

        ### Define the mode solver 
        ms = _2DSlab1LRhombusHole(h,Lz,b,e,vertices,
                                  num_bands,k_points,resolution,Mater,Envir)
        
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

if __name__ == '__name__':
    main()